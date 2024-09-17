use itertools::Itertools;
use minimizers::simd::packed::IntoBpIterator;

use crate::{MsSequence, Seq, Sketcher, SketcherBuilder};

/// A packed minimizer representation.
/// Bit width of the underlying alphabet is unspecified, and should not matter:
/// really this should only be used as a unique identifier of the kmer/minimizer
/// that is only checked for equality, not inspected.
type KmerVal = u64;
/// The start position of a minimizer.
type Pos = usize;

#[derive(Clone)]
pub struct MinimizerParams {
    pub k: usize,
    pub l: usize,
}

impl MinimizerParams {
    fn w(&self) -> usize {
        self.l - self.k + 1
    }

    fn minimizers<'s>(&'s self, seq: Seq<'s>) -> impl Iterator<Item = (Pos, KmerVal)> + 's {
        minimizers::simd::minimizer::minimizer_scalar_it::<false>(seq, self.k, self.w())
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
        let (min_poss, min_val): (Vec<Pos>, Vec<KmerVal>) = self.minimizers(seq).unzip();
        (
            MinimizerSketcher {
                params: self.clone(),
                min_poss,
            },
            // Simply reinterpret the u64 minimizer values as bytes on which the suffix array is built.
            MsSequence(
                min_val
                    .into_iter()
                    .flat_map(|x| x.to_ne_bytes())
                    .collect_vec(),
            ),
        )
    }
}

pub struct MinimizerSketcher {
    params: MinimizerParams,
    min_poss: Vec<Pos>,
}

impl Sketcher for MinimizerSketcher {
    fn sketch(&self, seq: &[u8]) -> Option<(MsSequence, usize)> {
        let (min_poss, min_val): (Vec<Pos>, Vec<KmerVal>) = self.params.minimizers(seq).unzip();
        Some((
            // Simply reinterpret the u64 minimizer values as bytes on which the suffix array is built.
            MsSequence(
                min_val
                    .into_iter()
                    .flat_map(|x| x.to_ne_bytes())
                    .collect_vec(),
            ),
            *min_poss.first()?,
        ))
    }

    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize> {
        const ALIGNMENT: usize = std::mem::size_of::<KmerVal>();
        if ms_pos % ALIGNMENT != 0 {
            return None;
        }
        Some(self.min_poss[ms_pos / ALIGNMENT])
    }
}
