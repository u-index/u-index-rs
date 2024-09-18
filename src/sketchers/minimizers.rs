use std::collections::HashMap;

use itertools::Itertools;
use mem_dbg::{MemDbg, MemSize};
use minimizers::simd::packed::IntoBpIterator;
use sux::traits::IndexedSeq;
use tracing::trace;

use crate::{utils::Timer, MsSequence, Seq, SketchError, Sketcher, SketcherBuilder};

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
}

impl MinimizerParams {
    fn w(&self) -> usize {
        self.l - self.k + 1
    }

    fn minimizers<'s>(&'s self, seq: Seq<'s>) -> impl Iterator<Item = (Pos, KmerVal)> + 's {
        minimizers::simd::minimizer::minimizer_simd_it::<false>(seq, self.k, self.w())
            .dedup()
            .map(move |pos| {
                let pos = pos as usize;
                (pos, (&&seq[pos..pos + self.k]).to_word() as KmerVal)
            })
    }
}

impl SketcherBuilder for MinimizerParams {
    type Sketcher = MinimizerSketcher;

    fn sketch(&self, seq: Seq) -> (Self::Sketcher, MsSequence) {
        assert!(
            self.k <= KmerVal::BITS as usize / 8,
            "k={} is too large to fit k bytes in a u64",
            self.k
        );
        trace!("Sequence length {}", seq.len());
        let mut timer = Timer::new("Computing minimizers");
        let (min_poss, min_val): (Vec<Pos>, Vec<KmerVal>) = self.minimizers(seq).unzip();
        trace!("Num minimizers: {}", min_poss.len());
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
            // When there is only a unique kmer as minimizer, use at least 1 byte still.
            let kmer_width = id.next_power_of_two().trailing_zeros().div_ceil(8).max(1) as usize;
            (kmer_map, kmer_width)
        } else {
            (HashMap::new(), self.k.min(8))
        };
        trace!("kmer_width: {kmer_width} bytes");
        timer.next("Building EF");
        let mut min_poss_ef = sux::dict::elias_fano::EliasFanoBuilder::new(
            min_poss.len(),
            *min_poss.last().unwrap_or(&0),
        );
        for p in min_poss {
            min_poss_ef.push(p);
        }
        let sketcher = MinimizerSketcher {
            params: self.clone(),
            min_poss: min_poss_ef.build_with_seq(),
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
                let buf = val.to_ne_bytes();
                bytes.extend_from_slice(&buf[..self.kmer_width]);
            }
        } else {
            for val in min_vals {
                let buf = val.to_ne_bytes();
                bytes.extend_from_slice(&buf[..self.kmer_width]);
            }
        }
        Some(MsSequence(bytes))
    }
}

impl Sketcher for MinimizerSketcher {
    fn sketch(&self, seq: Seq) -> Result<(MsSequence, usize), SketchError> {
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
        Some(self.min_poss.get(ms_pos / self.kmer_width))
    }
}
