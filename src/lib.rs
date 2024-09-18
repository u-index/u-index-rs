use std::cell::RefCell;

use indices::{IndexBuilderEnum, IndexEnum};
use mem_dbg::{MemDbg, MemSize};
use sketchers::{SketcherBuilderEnum, SketcherEnum};
use utils::{Timer, INIT_TRACE};

pub mod indices;
pub mod sketchers;
mod utils;

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
// - Currently we only take `u8` encoded 2-bit alphabet as input. That should be extended to
//   2-bit encoded data and full 8-bit alphabets.
// TODO: Generalize underlying KmerVal and Pos types.
// TODO: Streaming input
// TODO: Rolling hash for O(1) checking.
// TODO: True constant time checking?
// TODO: Consider whether `Index` should own the input.
// TODO: minmers?
// TODO: Randomize the minimizer order
// TODO: Perfect hash the minimizers
// TODO: Reserve ms ID for sequence ends.
// TODO: Use EF for sequence boundaries.
// TODO: Count false positives and other events. Add tracing or return them directly.

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

pub trait Index: MemSize + MemDbg {
    /// Return all places where the pattern occurs.
    fn query<'i>(&'i self, pattern: &[u8]) -> Box<dyn Iterator<Item = usize> + 'i>;
}

/// Sketch a plain sequence to minimizer space.
pub trait SketcherBuilder {
    type Sketcher: Sketcher + 'static;
    /// Take an input text, find its minimizers, and compress to the target space.
    fn sketch(&self, seq: Seq) -> (Self::Sketcher, MsSequence);
}

pub enum SketchError {
    /// The sequence is too short to be sketched, and can not be searched.
    TooShort,
    /// The sequence was sketched, but contains unknown minimizers, and thus is
    /// surely not contained in the text.
    UnknownMinimizer,
}

pub trait Sketcher: MemSize + MemDbg {
    /// Take an input text, compute its minimizers, and compress those into the
    /// target `u8` alphabet. This could be done a few ways, e.g.:
    /// - concatenating the KmerVals,
    /// - using a hash function to map the KmerVals to a smaller range.
    /// Also returns the position in `seq` of the first minimizer.
    /// Returns `None` when `seq` is too short to contain a minimizer.
    fn sketch(&self, seq: Seq) -> Result<(MsSequence, usize), SketchError>;

    /// Take a position of a character in the minimizer space, and return its start position in the original sequence.
    /// Returns `None` when the position in the minimizer space text is not aligned with the size of the encoded minimizers.
    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize>;
}

#[derive(MemSize, MemDbg)]
pub struct UIndex<'s> {
    seq: Seq<'s>,
    sketcher: SketcherEnum,
    ms_index: IndexEnum,
    query_stats: RefCell<QueryStats>,
}

#[derive(MemSize, MemDbg, Default, Debug)]
struct QueryStats {
    /// Pattern is too short to sketch.
    too_short: usize,
    /// Pattern contained a minimizer that is not in the input text.
    unknown_minimizer: usize,
    /// A minimizer-space match was found that does not align to the minimizer alphabet size.
    misaligned_ms_pos: usize,
    /// A minimizer-space match was found, but its containing sequence is too short to contain the pattern.
    out_of_bounds: usize,
    /// Mismatch in sequence space.
    mismatches: usize,
    /// Matches.
    matches: usize,
}

impl<'s> Drop for UIndex<'s> {
    fn drop(&mut self) {
        let stats = self.query_stats.borrow();
        tracing::info!("Query stats: {:#?}", *stats);
    }
}

impl<'s> UIndex<'s> {
    /// 1. Sketch input to minimizer space.
    /// 2. Build minimizer space index.
    pub fn build(
        seq: Seq<'s>,
        sketch_params: SketcherBuilderEnum,
        index_params: IndexBuilderEnum,
    ) -> Self {
        *INIT_TRACE;
        let mut timer = Timer::new("Sketch");
        let (sketcher, ms_seq) = sketch_params.sketch(seq);
        timer.next("Build");
        let ms_index = index_params.build(ms_seq.0);
        Self {
            seq,
            sketcher,
            ms_index,
            query_stats: RefCell::new(QueryStats::default()),
        }
    }

    /// 1. Sketch query to minimizer space.
    /// 2. Query the minimizer space index.
    /// 3. Check all occurrences.
    /// Returns `None` if the pattern is too short to contain a minimizer.
    /// When the pattern contains an unknown minimizer, an empty iterator is returned.
    pub fn query<'p>(&'p self, pattern: Seq<'p>) -> Option<Box<dyn Iterator<Item = usize> + 'p>> {
        let (ms_pattern, offset) = match self.sketcher.sketch(pattern) {
            Ok(x) => x,
            Err(SketchError::TooShort) => {
                self.query_stats.borrow_mut().too_short += 1;
                return None;
            }
            Err(SketchError::UnknownMinimizer) => {
                self.query_stats.borrow_mut().unknown_minimizer += 1;
                return Some(Box::new(std::iter::empty()));
            }
        };
        let ms_occ = self.ms_index.query(&ms_pattern.0);
        Some(Box::new(ms_occ.filter_map(move |ms_pos| {
            // Checking:
            // 1. Map minimizer space pos back to original space.
            // 2. Check query against original sequence.

            // The sketcher returns None when `ms_pos` does not align with a minimizer position.
            // The `checked_sub` fails when the minimizer is very close to the
            // start, and the `offset` doesn't fit before.
            let plain_pos = self.sketcher.ms_pos_to_plain_pos(ms_pos).or_else(|| {
                self.query_stats.borrow_mut().misaligned_ms_pos += 1;
                None
            })?;
            let pos = plain_pos.checked_sub(offset).or_else(|| {
                self.query_stats.borrow_mut().out_of_bounds += 1;
                None
            })?;

            // The `get` fails when the minimizer match is too close to the end
            // and the pattern doesn't fit after it.
            let matches = self.seq.get(pos..pos + pattern.len()).or_else(|| {
                self.query_stats.borrow_mut().out_of_bounds += 1;
                None
            })? == pattern;
            if matches {
                self.query_stats.borrow_mut().matches += 1;
                Some(pos)
            } else {
                self.query_stats.borrow_mut().mismatches += 1;
                None
            }
        })))
    }
}

#[cfg(test)]
mod test {
    use indices::DivSufSortSa;
    use mem_dbg::SizeFlags;
    use sketchers::{IdentityParams, MinimizerParams};
    use tracing::trace;

    use super::*;

    #[test]
    fn test_identity_simple() {
        let seq = b"ACGTACGTACGTACGT";
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa);
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
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa);
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
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa);
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

        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa);
        let index = UIndex::build(&seq, sketcher, ms_index);

        for remap in [false, true] {
            for l in [1, 10, 100] {
                for k in [1, 2, 3, 4, 5, 6, 7, 8] {
                    if k > l {
                        continue;
                    }
                    let sketcher =
                        SketcherBuilderEnum::Minimizer(sketchers::MinimizerParams { l, k, remap });
                    let uindex = UIndex::build(&seq, sketcher, ms_index);
                    for _ in 0..1000 {
                        let len = l + rand::random::<usize>() % 100;
                        let pos = rand::random::<usize>() % (seq.len() - len);
                        let query = &seq[pos..pos + len];

                        let mut index_occ = index.query(query).unwrap().collect::<Vec<_>>();
                        let mut uindex_occ = uindex.query(query).unwrap().collect::<Vec<_>>();
                        index_occ.sort();
                        uindex_occ.sort();
                        assert_eq!(
                            index_occ, uindex_occ,
                            "l {l} k {k} remap {remap} pos {pos} query {query:?}"
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn test_minspace_negative() {
        let seq = (0..1000000)
            .map(|_i| rand::random::<u8>() % 4)
            .collect::<Vec<_>>();

        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa);
        let index = UIndex::build(&seq, sketcher, ms_index);

        for remap in [false, true] {
            for l in [1, 10, 100] {
                for k in [1, 2, 3, 4, 5, 6, 7, 8] {
                    if k > l {
                        continue;
                    }
                    let sketcher =
                        SketcherBuilderEnum::Minimizer(sketchers::MinimizerParams { l, k, remap });
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

    fn read_human_genome() -> Sequence {
        let mut timer = Timer::new("Reading");
        let Ok(mut reader) = needletail::parse_fastx_file("human-genome.fa") else {
            panic!("Did not find human-genome.fa. Add/symlink it to test runtime on it.");
        };
        let mut seq = vec![];
        while let Some(r) = reader.next() {
            let r = r.unwrap();
            seq.extend_from_slice(&r.seq());
            break;
        }
        timer.next("Mapping");
        for c in &mut seq {
            *c = match *c {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => panic!("Unexpected character {c}"),
            }
        }
        seq
    }

    #[test]
    #[ignore = "needs human-genome.fa"]
    fn human_genome() {
        let seq = read_human_genome();

        for l in [200] {
            for k in [5, 6, 7, 8] {
                for remap in [true] {
                    if k > l {
                        continue;
                    }

                    let mut timer = Timer::new("Build");
                    trace!("remap {remap} l {l} k {k}");
                    let sketcher = SketcherBuilderEnum::Minimizer(MinimizerParams { l, k, remap });
                    let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa);
                    let uindex = UIndex::build(&seq, sketcher, ms_index);
                    trace!(
                        "Index size:   {} MB",
                        uindex.mem_size(SizeFlags::default()) / (1024 * 1024)
                    );
                    timer.next("Query");
                    for _ in 0..1000 {
                        let len = 2 * l + rand::random::<usize>() % l;
                        let pos = rand::random::<usize>() % (seq.len() - len);
                        let query = &seq[pos..pos + len];

                        let uindex_occ = uindex.query(query).unwrap().collect::<Vec<_>>();
                        assert!(uindex_occ.contains(&pos));
                    }
                }
            }
        }
    }
}
