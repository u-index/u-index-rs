use std::collections::HashMap;

use cacheline_ef::CachelineEfVec;
use itertools::Itertools;
use mem_dbg::{MemDbg, MemSize};
use minimizers::Captures;
use packed_seq::Seq;
use sux::traits::IndexedSeq;
use tracing::trace;

use crate::{
    utils::{Stats, Timer},
    MsSequence, SketchError, Sketcher, SketcherBuilder, S,
};

/// A packed minimizer representation.
/// Bit width of the underlying alphabet is unspecified, and should not matter:
/// really this should only be used as a unique identifier of the kmer/minimizer
/// that is only checked for equality, not inspected.
type KmerVal = u64;
/// The start position of a minimizer.
type Pos = usize;

#[derive(Clone, Copy, MemSize, MemDbg)]
pub struct MinimizerParams {
    /// The kmer/minimizer size.
    pub k: usize,
    /// Each sequence of length l will contain at least one minimizer.
    pub l: usize,
    /// When true, remap kmers to a smaller alphabet.
    pub remap: bool,
    /// Which implementation to use for the minimizer positions.
    pub cacheline_ef: bool,
}

impl MinimizerParams {
    fn w(&self) -> usize {
        self.l - self.k + 1
    }

    fn minimizers<'s>(
        &'s self,
        seq: impl Seq<'s>,
    ) -> impl Iterator<Item = (Pos, KmerVal)> + Captures<&'s ()> {
        minimizers::simd::minimizer::minimizer_scalar_it::<false>(seq, self.k, self.w())
            .dedup()
            .map(move |pos| {
                let pos = pos as usize;
                (
                    pos,
                    packed_seq::Seq::to_word(&seq.slice(pos..pos + self.k)) as KmerVal,
                )
            })
    }

    fn minimizers_par<'s>(
        &'s self,
        seq: impl Seq<'s>,
    ) -> impl Iterator<Item = (Pos, KmerVal)> + Captures<&'s ()> {
        minimizers::simd::minimizer::minimizer_simd_it::<false>(seq, self.k, self.w())
            .dedup()
            .map(move |pos| {
                let pos = pos as usize;
                (
                    pos,
                    packed_seq::Seq::to_word(&seq.slice(pos..pos + self.k)) as KmerVal,
                )
            })
    }
}

impl SketcherBuilder for MinimizerParams {
    type Sketcher = MinimizerSketcher;

    fn sketch_with_stats(&self, seq: S<'_>, stats: &Stats) -> (Self::Sketcher, MsSequence) {
        assert!(
            self.k <= KmerVal::BITS as usize / 2,
            "k={} is too large to fit k bytes in a u64",
            self.k
        );
        trace!("Sequence length {}", seq.len());
        stats.set("sequence_length", seq.len());
        let mut timer = Timer::new_stats("computing_minimizers", stats);
        let (min_poss, min_val): (Vec<Pos>, Vec<KmerVal>) = self.minimizers_par(seq).unzip();
        trace!("Num minimizers: {}", min_poss.len());
        stats.set("num_minimizers", min_poss.len());
        let (kmer_map, kmer_width) = if self.remap {
            timer.next("Building remap");
            let mut kmer_map = HashMap::new();
            let mut id = 0usize;
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
            stats.set("kmer_width_bits", 2 * self.k);
            (HashMap::new(), self.k.div_ceil(4))
        };
        trace!("kmer_width: {kmer_width} bytes");
        stats.set("kmer_width", kmer_width);
        timer.next("Building EF");
        let min_poss_ef = if self.cacheline_ef {
            let min_poss_ef = sux::dict::elias_fano::EliasFanoBuilder::new(0, 0);
            min_poss_ef.build_with_seq()
        } else {
            let mut min_poss_ef = sux::dict::elias_fano::EliasFanoBuilder::new(
                min_poss.len(),
                *min_poss.last().unwrap_or(&0),
            );
            for &p in &min_poss {
                min_poss_ef.push(p);
            }
            min_poss_ef.build_with_seq()
        };
        let min_poss_cacheline_ef = if self.cacheline_ef {
            CachelineEfVec::new(&min_poss.into_iter().map(|x| x as u64).collect_vec())
        } else {
            CachelineEfVec::new(&vec![])
        };
        let sketcher = MinimizerSketcher {
            params: self.clone(),
            min_poss: min_poss_ef,
            min_poss_cacheline_ef,
            kmer_map,
            kmer_width,
        };
        timer.next("Remapping");
        let ms_sequence = sketcher
            .remap_minimizer_values(&min_val)
            .expect("All minimizers of the input should be found");
        (sketcher, ms_sequence)
    }
}

#[derive(MemSize, MemDbg)]
pub struct MinimizerSketcher {
    params: MinimizerParams,
    /// Positions in the plain sequence of all minimizers.
    min_poss: sux::dict::elias_fano::EfSeq,
    /// Positions in the plain sequence of all minimizers.
    min_poss_cacheline_ef: cacheline_ef::CachelineEfVec,
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

impl Sketcher for MinimizerSketcher {
    fn width(&self) -> usize {
        self.kmer_width
    }

    fn k(&self) -> usize {
        self.params.k
    }

    /// Return the number of minimizers.
    fn len(&self) -> usize {
        self.min_poss.len()
    }

    fn sketch(&self, seq: S) -> Result<(MsSequence, usize), SketchError> {
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
        if self.params.cacheline_ef {
            Some(self.min_poss_cacheline_ef.index(ms_pos / self.kmer_width) as usize)
        } else {
            Some(self.min_poss.get(ms_pos / self.kmer_width))
        }
    }

    fn get_ms_minimizer(&self, ms_seq: &[u8], ms_pos: usize) -> Option<usize> {
        let mut val = [0u8; 8];
        val[8 - self.kmer_width..].copy_from_slice(&ms_seq[ms_pos..ms_pos + self.kmer_width]);
        Some(usize::from_be_bytes(val))
    }

    fn get_ms_minimizer_via_plaintext(&self, seq: S, ms_pos: usize) -> Option<usize> {
        let pos = self.ms_pos_to_plain_pos(ms_pos)?;
        let kmer = packed_seq::Seq::to_word(&seq.slice(pos..pos + self.params.k)) as KmerVal;
        if self.params.remap {
            self.kmer_map.get(&kmer).copied()
        } else {
            Some(kmer as usize)
        }
    }
}
