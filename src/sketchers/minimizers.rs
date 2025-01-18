use std::collections::HashMap;

use cacheline_ef::CachelineEfVec;
use itertools::Itertools;
use mem_dbg::MemSize;
use packed_seq::{Seq, SeqVec};
use serde_json::Value;
use sux::traits::IndexedSeq;
use tracing::{info, trace};

use crate::{
    utils::{Stats, Timer},
    MsSequence, SketchError, Sketcher, SketcherBuilder,
};

/// A packed minimizer representation.
/// Bit width of the underlying alphabet is unspecified, and should not matter:
/// really this should only be used as a unique identifier of the kmer/minimizer
/// that is only checked for equality, not inspected.
type KmerVal = u64;
/// The start position of a minimizer.
type Pos = usize;

#[derive(Clone, Copy, MemSize, Debug)]
pub struct MinimizerParams {
    /// The kmer/minimizer size.
    pub k: usize,
    /// Each sequence of length l will contain at least one minimizer.
    pub l: usize,
    /// When true, remap kmers to a smaller alphabet.
    pub remap: bool,
    /// Which implementation to use for the minimizer positions.
    pub cacheline_ef: bool,
    /// When true, remapping skips 0.
    pub skip_zero: bool,
}

impl MinimizerParams {
    fn w(&self) -> usize {
        self.l - self.k + 1
    }

    fn minimizers<'s, S: Seq<'s>>(
        &self,
        seq: S,
    ) -> impl Iterator<Item = (Pos, KmerVal)> + use<'s, S> {
        let k = self.k;
        let w = self.w();

        // FIXME: Use SIMD implementation?
        // FIXME: Use iterator version?
        let mut out = vec![];
        simd_minimizers::minimizer_positions_scalar(seq, k, w, &mut out);
        out.into_iter().map(move |pos| {
            let pos = pos as usize;
            (
                pos,
                packed_seq::Seq::to_word(&seq.slice(pos..pos + k)) as KmerVal,
            )
        })
        // let mut out = vec![];
        // minimizers::simd::minimizer::minimizers_collect_and_dedup::<false>(seq, k, w, &mut out);
        // out.into_iter().map(move |pos| {
        //     let pos = pos as usize;
        //     (
        //         pos,
        //         packed_seq::Seq::to_word(&seq.slice(pos..pos + k)) as KmerVal,
        //     )
        // })
    }

    fn minimizers_par<'s, S: Seq<'s>>(
        &self,
        seq: S,
    ) -> impl Iterator<Item = (Pos, KmerVal)> + use<'s, S> {
        let k = self.k;
        let w = self.w();
        let mut out = vec![];
        simd_minimizers::minimizer_positions(seq, k, w, &mut out);
        out.into_iter().map(move |pos| {
            let pos = pos as usize;
            (
                pos,
                packed_seq::Seq::to_word(&seq.slice(pos..pos + k)) as KmerVal,
            )
        })
    }
}

impl<SV: SeqVec> SketcherBuilder<SV> for MinimizerParams {
    fn sketch_with_stats(
        &self,
        seq: SV::Seq<'_>,
        stats: &Stats,
    ) -> (Box<dyn Sketcher<SV>>, MsSequence) {
        stats.set_val("sketcher", Value::String("Minimizers".to_string()));
        stats.set("sketch_k", self.k);
        stats.set("sketch_l", self.l);
        stats.set("sketch_remap", self.remap as u64);
        stats.set("sketch_cacheline_ef", self.cacheline_ef as u64);
        stats.set("sketch_skip_zero", self.skip_zero as u64);

        assert!(
            self.k <= KmerVal::BITS as usize / 2,
            "k={} is too large to fit k bytes in a u64",
            self.k
        );
        trace!("Sequence length {}", seq.len());
        trace!("k: {}", self.k);
        trace!("l: {}", self.l);
        trace!("remap: {}", self.remap);
        trace!("cacheline_ef: {}", self.cacheline_ef);
        stats.set("sequence_length", seq.len());
        let mut timer = Timer::new_stats("computing_minimizers", stats);
        let (min_poss, min_val): (Vec<Pos>, Vec<KmerVal>) = self.minimizers_par(seq).unzip();
        info!("Num minimizers: {:>9}", min_poss.len());
        stats.set("num_minimizers", min_poss.len());
        let (kmer_map, kmer_width) = if self.remap {
            timer.next("Building remap");
            let mut kmer_map = HashMap::new();
            // Start at 1 to avoid 0 bytes.
            let mut id = if self.skip_zero { 1usize } else { 0usize };
            for &kmer in &min_val {
                if !kmer_map.contains_key(&kmer) {
                    kmer_map.insert(kmer, id);
                    id += 1;
                }
            }
            trace!("Num distinct minimizers: {}", id);
            stats.set("num_distinct_minimizers", id);
            // When there is only a unique kmer as minimizer, use at least 1 byte still.
            let kmer_width_bits = id.next_power_of_two().trailing_zeros();
            trace!("kmer_width: {kmer_width_bits} bits");
            stats.set("kmer_width_bits", kmer_width_bits);
            let kmer_width = kmer_width_bits.div_ceil(8).max(1) as usize;
            (kmer_map, kmer_width)
        } else {
            trace!(
                "Using raw minimizer values. Max: {}",
                min_val.iter().max().unwrap_or(&0)
            );
            stats.set("kmer_width_bits", 2 * self.k);
            (HashMap::new(), self.k.div_ceil(4))
        };
        info!("kmer_width: {kmer_width} bytes");
        stats.set("kmer_width", kmer_width);
        timer.next("Building EF");
        let min_poss = if self.cacheline_ef {
            MinimizerPositions::CachelineEf(CachelineEfVec::new(
                &min_poss.into_iter().map(|x| x as u64).collect_vec(),
            ))
        } else {
            let mut min_poss_ef = sux::dict::elias_fano::EliasFanoBuilder::new(
                min_poss.len(),
                *min_poss.last().unwrap_or(&0),
            );
            for &p in &min_poss {
                min_poss_ef.push(p);
            }
            MinimizerPositions::EliasFano(min_poss_ef.build_with_seq())
        };
        let sketcher = MinimizerSketcher {
            params: self.clone(),
            min_poss,
            kmer_map,
            kmer_width,
        };
        timer.next("Remapping");
        let ms_sequence = sketcher
            .remap_minimizer_values(&min_val)
            .expect("All minimizers of the input should be found");
        (Box::new(sketcher), ms_sequence)
    }
}

#[derive(MemSize)]
enum MinimizerPositions {
    EliasFano(sux::dict::elias_fano::EfSeq),
    CachelineEf(cacheline_ef::CachelineEfVec),
}

#[derive(MemSize)]
pub struct MinimizerSketcher {
    params: MinimizerParams,
    /// Positions in the plain sequence of all minimizers.
    min_poss: MinimizerPositions,
    /// When `remap` is true, a map from kmers to smaller IDs.
    kmer_map: HashMap<KmerVal, usize>,
    /// The width in bytes of the kmer IDs used.
    kmer_width: usize,
}

impl MinimizerSketcher {
    /// Takes minimizer values, and returns their byte representation.
    /// Returns `None` when a minimizer was found that is not a minimizer in the input.
    fn remap_minimizer_values(&self, min_vals: &[KmerVal]) -> Option<MsSequence> {
        let mut bytes = Vec::with_capacity(min_vals.len() * self.kmer_width);
        if self.params.remap {
            for val in min_vals {
                let val = *self.kmer_map.get(val)?;
                let buf = val.to_be_bytes();
                bytes.extend_from_slice(&buf[8 - self.kmer_width..]);
            }
        } else {
            for val in min_vals {
                let buf = val.to_be_bytes();
                bytes.extend_from_slice(&buf[8 - self.kmer_width..]);
            }
        }
        Some(MsSequence(bytes))
    }
}

impl<SV: SeqVec> Sketcher<SV> for MinimizerSketcher {
    fn width(&self) -> usize {
        self.kmer_width
    }

    fn k(&self) -> usize {
        self.params.k
    }

    /// Return the number of minimizers.
    fn len(&self) -> usize {
        match &self.min_poss {
            MinimizerPositions::EliasFano(elias_fano) => elias_fano.len(),
            MinimizerPositions::CachelineEf(cacheline_ef_vec) => cacheline_ef_vec.len(),
        }
    }

    fn sketch(&self, seq: SV::Seq<'_>) -> Result<(MsSequence, usize), SketchError> {
        let (min_poss, min_vals): (Vec<Pos>, Vec<KmerVal>) = self.params.minimizers(seq).unzip();
        let offset = *min_poss.first().ok_or(SketchError::TooShort)?;
        Ok((
            self.remap_minimizer_values(&min_vals)
                .ok_or(SketchError::UnknownMinimizer)?,
            offset,
        ))
    }

    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize> {
        if ms_pos % self.kmer_width != 0 {
            return None;
        }
        match &self.min_poss {
            MinimizerPositions::EliasFano(elias_fano) => {
                Some(elias_fano.get(ms_pos / self.kmer_width))
            }
            MinimizerPositions::CachelineEf(cacheline_ef_vec) => {
                Some(cacheline_ef_vec.index(ms_pos / self.kmer_width) as usize)
            }
        }
    }

    fn get_ms_minimizer(&self, ms_seq: &[u8], ms_pos: usize) -> Option<usize> {
        let mut val = [0u8; 8];
        val[8 - self.kmer_width..].copy_from_slice(&ms_seq[ms_pos..ms_pos + self.kmer_width]);
        Some(usize::from_be_bytes(val))
    }

    fn get_ms_minimizer_via_plaintext(&self, seq: SV::Seq<'_>, ms_pos: usize) -> Option<usize> {
        let pos = Sketcher::<SV>::ms_pos_to_plain_pos(self, ms_pos)?;
        let kmer = seq.slice(pos..pos + self.params.k).to_word() as KmerVal;
        if self.params.remap {
            self.kmer_map.get(&kmer).copied()
        } else {
            Some(kmer as usize)
        }
    }
}
