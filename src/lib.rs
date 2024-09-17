pub mod indices;
pub mod sketchers;

// Terminology and variables:
// - Index: a datastructure that returns all locations where a pattern matches.
// - k-mer and minimizer are the used interchangeably.
// - k: the k-mer length of the minimizers.
// - l: the length of guaranteed preserved matches.
// - w = l - k + 1: the number of k-mers in a preserved window.
// - seq: the input sequence.
//   NOTE: Avoid 'text', always use 'seq'.
// - MS: MinSpace
// - ms_seq: the minimizer-space compressed text.
// - plain seq: the original sequence.

// TODO: Generalize the alphabet.
// - Currently the underlying index is assumed to only support `u8` alphabet.
// - Currently we only take `u8` alphabet as input. That should be extended to
//   2-bit encoded data.
// TODO: Generalize underlying KmerVal and Pos types.
// TODO: Streaming input
// TODO: Rolling hash for O(1) checking.
// TODO: True constant time checking?
// TODO: Consider whether `Index` should own the input.
// TODO: Compress the minimizer alphabet.
// TODO: Use EF for minimizer locations.

/// A reference to a sequence over a u8 alphabet.
type Seq<'s> = &'s [u8];
/// A sequence over a u8 alphabet.
type Sequence = Vec<u8>;
/// A reference to a minimizer space sequence.
// struct MsSeq<'s>(&'s [u8]);
/// A minimizer space sequence.
pub struct MsSequence(Vec<u8>);

/// A generic index to locate strings.
/// The index owns the input text.
pub trait IndexBuilder {
    type Index: Index + 'static;
    /// Build an index on the text.
    fn build(&self, text: Sequence) -> Self::Index;
}

pub trait Index {
    /// Return all places where the pattern occurs.
    fn query<'i>(&'i self, pattern: &[u8]) -> Box<dyn Iterator<Item = usize> + 'i>;
}

/// Sketch a plain sequence to minimizer space.
pub trait SketcherBuilder {
    type Sketcher: Sketcher + 'static;
    /// Take an input text, find its minimizers, and compress to the target space.
    fn sketch(&self, seq: Seq) -> (Self::Sketcher, MsSequence);
}

pub trait Sketcher {
    /// Take an input text, compute its minimizers, and compress those into the
    /// target `u8` alphabet. This could be done a few ways, e.g.:
    /// - concatenating the KmerVals,
    /// - using a hash function to map the KmerVals to a smaller range.
    /// Also returns the position in `seq` of the first minimizer.
    /// Returns `None` when `seq` is too short to contain a minimizer.
    fn sketch(&self, seq: Seq) -> Option<(MsSequence, usize)>;

    /// Take a position of a character in the minimizer space, and return its start position in the original sequence.
    /// Returns `None` when the position in the minimizer space text is not aligned with the size of the encoded minimizers.
    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize>;
}

pub struct UIndex<'s> {
    seq: Seq<'s>,
    sketcher: Box<dyn Sketcher>,
    ms_index: Box<dyn Index>,
}

impl<'s> UIndex<'s> {
    /// 1. Sketch input to minimizer space.
    /// 2. Build minimiser space sequence.
    /// 2. (optionally) remap to smaller alphabet
    /// 3. Concat the minimizer values.
    /// 4. Reinterpret as underlying 2bit alphabet.
    /// 5. Build underlying index.
    pub fn build(
        seq: Seq<'s>,
        sketch_params: impl SketcherBuilder,
        index_params: impl IndexBuilder,
    ) -> Self {
        let (sketcher, ms_seq) = sketch_params.sketch(seq);
        let ms_index = Box::new(index_params.build(ms_seq.0));
        Self {
            seq,
            sketcher: Box::new(sketcher),
            ms_index,
        }
    }

    /// 1. Sketch query to minimizer space.
    /// 2. Query the minimizer space index.
    /// 3. Check all occurrences.
    /// Returns `None` if the pattern is too short and does not contain a minimizer.
    pub fn query<'p>(&'p self, pattern: Seq<'p>) -> Option<Box<dyn Iterator<Item = usize> + 'p>> {
        let (ms_pattern, offset) = self.sketcher.sketch(pattern)?;
        let ms_occ = self.ms_index.query(&ms_pattern.0);
        Some(Box::new(ms_occ.filter_map(move |ms_pos| {
            // Checking:
            // 1. Map minimizer space pos back to original space.
            // 2. Check query against original sequence.

            // The sketcher returns None when `ms_pos` does not align with a minimizer position.
            // The `checked_sub` fails when the minimizer is very close to the
            // start, and the `offset` doesn't fit before.
            let pos = self
                .sketcher
                .ms_pos_to_plain_pos(ms_pos)?
                .checked_sub(offset)?;

            // The `get` fails when the minimizer match is too close to the end
            // and the pattern doesn't fit after it.
            let matches = self.seq.get(pos..pos + pattern.len())? == pattern;
            matches.then_some(pos)
        })))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_identity_simple() {
        let seq = b"ACGTACGTACGTACGT";
        let sketcher = sketchers::Identity;
        let ms_index = indices::DivSufSortSa;
        let uindex = UIndex::build(seq, sketcher, ms_index);
        let query = b"ACGT";
        let mut occ = uindex.query(query).unwrap().collect::<Vec<_>>();
        occ.sort();
        assert_eq!(occ, vec![0, 4, 8, 12]);
    }
    #[test]
    fn test_identity_positive() {
        let seq = (0..1000000)
            .map(|_i| rand::random::<u8>())
            .collect::<Vec<_>>();
        let sketcher = sketchers::Identity;
        let ms_index = indices::DivSufSortSa;
        let uindex = UIndex::build(&seq, sketcher, ms_index);
        for _ in 0..1000 {
            let len = rand::random::<usize>() % 100;
            let pos = rand::random::<usize>() % (seq.len() - len);
            let query = &seq[pos..pos + len];
            let occ = uindex.query(query).unwrap().collect::<Vec<_>>();
            assert!(occ.len() > 0);
            for &pos in &occ {
                assert_eq!(&seq[pos..pos + len], query);
            }
        }
    }
    #[test]
    fn test_identity_negative() {
        let seq = (0..1000000)
            .map(|_i| rand::random::<u8>())
            .collect::<Vec<_>>();
        let sketcher = sketchers::Identity;
        let ms_index = indices::DivSufSortSa;
        let uindex = UIndex::build(&seq, sketcher, ms_index);
        for _ in 0..1000 {
            let len = 6;
            let query = (0..len).map(|_i| rand::random::<u8>()).collect::<Vec<_>>();
            let occ = uindex.query(&query).unwrap().collect::<Vec<_>>();
            assert_eq!(occ.len(), 0);
        }
    }
    #[test]
    fn test_minspace_positive() {
        let seq = (0..1000000)
            .map(|_i| rand::random::<u8>() % 4)
            .collect::<Vec<_>>();

        let ms_index = indices::DivSufSortSa;
        let index = UIndex::build(&seq, sketchers::Identity, ms_index);
        for l in [1, 10, 100] {
            for k in [1, 2, 3, 4, 5, 6, 7, 8] {
                if k > l {
                    continue;
                }
                let sketcher = sketchers::MinimizerParams { l, k };
                let uindex = UIndex::build(&seq, sketcher, ms_index);
                for _ in 0..1000 {
                    let len = l + rand::random::<usize>() % 100;
                    let pos = rand::random::<usize>() % (seq.len() - len);
                    let query = &seq[pos..pos + len];

                    let mut index_occ = index.query(query).unwrap().collect::<Vec<_>>();
                    let mut uindex_occ = uindex.query(query).unwrap().collect::<Vec<_>>();
                    index_occ.sort();
                    uindex_occ.sort();
                    assert_eq!(index_occ, uindex_occ);
                }
            }
        }
    }
    #[test]
    fn test_minspace_negative() {
        let seq = (0..1000000)
            .map(|_i| rand::random::<u8>() % 4)
            .collect::<Vec<_>>();

        let ms_index = indices::DivSufSortSa;
        let index = UIndex::build(&seq, sketchers::Identity, ms_index);

        for l in [1, 10, 100] {
            for k in [1, 2, 3, 4, 5, 6, 7, 8] {
                if k > l {
                    continue;
                }
                let sketcher = sketchers::MinimizerParams { l, k };
                let uindex = UIndex::build(&seq, sketcher, ms_index);
                for _ in 0..1000 {
                    let len = l + rand::random::<usize>() % 100;
                    let query = &(0..len)
                        .map(|_i| rand::random::<u8>() % 4)
                        .collect::<Vec<_>>();

                    let mut index_occ = index.query(query).unwrap().collect::<Vec<_>>();
                    let mut uindex_occ = uindex.query(query).unwrap().collect::<Vec<_>>();
                    index_occ.sort();
                    uindex_occ.sort();
                    assert_eq!(index_occ, uindex_occ);
                }
            }
        }
    }
}
