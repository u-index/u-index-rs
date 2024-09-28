use std::{cell::RefCell, collections::HashMap};

use indices::{IndexBuilderEnum, IndexEnum};
use mem_dbg::{MemDbg, MemSize, SizeFlags};
use packed_seq::{AsciiSeq, AsciiSeqVec, Seq, SeqVec};
use pyo3::pyclass;
use sketchers::{SketcherBuilderEnum, SketcherEnum};
use sux::traits::SuccUnchecked;
use tracing::trace;
use utils::{Stats, Timer, INIT_TRACE};

pub mod indices;
pub mod sketchers;
pub mod utils;

mod bench;
#[cfg(feature = "python-bindings")]
mod py;

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

/// A minimizer space sequence.
pub struct MsSequence(Vec<u8>);

/// A generic index to locate strings.
/// The index owns the input text.
pub trait IndexBuilder {
    type Index: Index + 'static;

    /// Build an index on the text, and keep track of statistics.
    /// `width` gives the width of each minimizer in bytes.
    /// Effectively, only a suffix array on the `width` byte wide integers is needed.
    fn build_with_stats(&self, text: Vec<u8>, width: usize, stats: &Stats) -> Self::Index;

    /// Build an index on the text.
    fn build(&self, text: Vec<u8>, width: usize) -> Self::Index {
        self.build_with_stats(text, width, &Stats::default())
    }
}

pub trait Index: MemSize + MemDbg {
    /// Return all places where the pattern occurs.
    fn query<'i>(&'i self, pattern: &[u8]) -> Box<dyn Iterator<Item = usize> + 'i>;
}

/// Sketch a plain sequence to minimizer space.
pub trait SketcherBuilder {
    type Sketcher: Sketcher + 'static;
    /// Take an input text, find its minimizers, and compress to the target space.
    /// Additionally log statistics to `stats`.
    fn sketch_with_stats(&self, seq: AsciiSeq, stats: &Stats) -> (Self::Sketcher, MsSequence);

    /// Take an input text, find its minimizers, and compress to the target space.
    fn sketch(&self, seq: AsciiSeq) -> (Self::Sketcher, MsSequence) {
        self.sketch_with_stats(seq, &Stats::default())
    }
}

pub enum SketchError {
    /// The sequence is too short to be sketched, and can not be searched.
    TooShort,
    /// The sequence was sketched, but contains unknown minimizers, and thus is
    /// surely not contained in the text.
    UnknownMinimizer,
}

pub trait Sketcher: MemSize + MemDbg {
    /// Returns the width in bytes of each minimizer.
    fn width(&self) -> usize;

    /// Take an input text, compute its minimizers, and compress those into the
    /// target `u8` alphabet. This could be done a few ways, e.g.:
    /// - concatenating the KmerVals,
    /// - using a hash function to map the KmerVals to a smaller range.
    /// Also returns the position in `seq` of the first minimizer.
    /// Returns `None` when `seq` is too short to contain a minimizer.
    fn sketch(&self, seq: AsciiSeq) -> Result<(MsSequence, usize), SketchError>;

    /// Take a position of a character in the minimizer space, and return its start position in the original sequence.
    /// Returns `None` when the position in the minimizer space text is not aligned with the size of the encoded minimizers.
    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize>;
}

#[derive(MemSize, MemDbg)]
#[pyclass]
pub struct UIndex {
    seq: AsciiSeqVec,
    sketcher: SketcherEnum,
    ms_index: IndexEnum,
    query_stats: RefCell<QueryStats>,
    stats: Stats,
    ranges: sux::dict::elias_fano::EfDict,
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
    /// Bad ranges.
    bad_ranges: usize,

    /// Total time in ns of sketching queries.
    t_sketch: usize,
    /// Total time in ns of looking up queries in the minimizer space index.
    t_search: usize,
    /// Total time in ns of converting the ms position back to the plain position.
    t_invert_pos: usize,
    /// Total time in ns checking matches.
    t_check: usize,
    /// Total time in ns checking ranges.
    t_ranges: usize,
}

impl Drop for UIndex {
    fn drop(&mut self) {
        let stats = self.query_stats.borrow();
        tracing::info!("Query stats: {:#?}", *stats);
    }
}

impl UIndex {
    /// 1. Sketch input to minimizer space.
    /// 2. Build minimizer space index.
    pub fn build(
        mut seq: AsciiSeqVec,
        sketch_params: SketcherBuilderEnum,
        index_params: IndexBuilderEnum,
    ) -> Self {
        *INIT_TRACE;
        let stats = Stats::default();
        let mut timer = Timer::new_stats("Sketch", &stats);
        let (sketcher, ms_seq) = sketch_params.sketch_with_stats(seq.as_slice(), &stats);
        timer.next("Build");
        let ms_index = index_params.build_with_stats(ms_seq.0, sketcher.width(), &stats);
        drop(timer);

        // Build seq ranges.
        let mut ef_ranges = sux::dict::elias_fano::EliasFanoBuilder::new(
            2 * seq.ranges.len(),
            seq.ranges.last().unwrap().1,
        );
        for (start, end) in std::mem::take(&mut seq.ranges) {
            ef_ranges.push(start);
            ef_ranges.push(end);
        }

        let uindex = Self {
            seq,
            sketcher,
            ms_index,
            query_stats: RefCell::new(QueryStats::default()),
            stats,
            ranges: ef_ranges.build_with_dict(),
        };
        uindex.stats.add(
            "seq_size_MB",
            uindex.seq.mem_size(SizeFlags::default()) as f32 / 1000000.,
        );
        let sketch_size = uindex.sketcher.mem_size(SizeFlags::default()) as f32 / 1000000.;
        uindex.stats.add("sketch_size_MB", sketch_size);
        trace!("Sketch size:   {sketch_size:>8.3} MB",);
        let index_size = uindex.ms_index.mem_size(SizeFlags::default()) as f32 / 1000000.;
        uindex.stats.add("index_size_MB", index_size);
        trace!("Index  size:   {index_size:>8.3} MB",);
        let ranges_size = uindex.ranges.mem_size(SizeFlags::default()) as f32 / 1000000.;
        uindex.stats.add("ranges_size_MB", index_size);
        trace!("Ranges size:   {ranges_size:>8.3} MB",);
        let total_size = sketch_size + index_size + ranges_size;
        uindex.stats.add("total_size_MB", total_size);
        trace!("Total  size:   {total_size:>8.3} MB",);
        uindex
    }

    pub fn stats(&self) -> HashMap<&'static str, f32> {
        let stats = self.stats.clone();
        let qs = self.query_stats.borrow();
        stats.set("query_too_short", qs.too_short);
        stats.set("query_unknown_minimizer", qs.unknown_minimizer);
        stats.set("query_misaligned_ms_pos", qs.misaligned_ms_pos);
        stats.set("query_out_of_bounds", qs.out_of_bounds);
        stats.set("query_mismatches", qs.mismatches);
        stats.set("query_matches", qs.matches);
        stats.set("t_query_sketch", qs.t_sketch as f32 / 1_000_000_000.);
        stats.set("t_query_search", qs.t_search as f32 / 1_000_000_000.);
        stats.set(
            "t_query_invert_pos",
            qs.t_invert_pos as f32 / 1_000_000_000.,
        );
        stats.set("t_query_check", qs.t_check as f32 / 1_000_000_000.);

        stats.into()
    }

    /// 1. Sketch query to minimizer space.
    /// 2. Query the minimizer space index.
    /// 3. Check all occurrences.
    /// Returns `None` if the pattern is too short to contain a minimizer.
    /// When the pattern contains an unknown minimizer, an empty iterator is returned.
    pub fn query<'p>(
        &'p self,
        pattern: AsciiSeq<'p>,
    ) -> Option<Box<dyn Iterator<Item = usize> + 'p>> {
        let t1 = std::time::Instant::now();
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
        let t2 = std::time::Instant::now();
        self.query_stats.borrow_mut().t_sketch += t2.duration_since(t1).subsec_nanos() as usize;
        let ms_occ = self.ms_index.query(&ms_pattern.0);
        let t3 = std::time::Instant::now();
        self.query_stats.borrow_mut().t_search += t3.duration_since(t2).subsec_nanos() as usize;
        let mut last = t3;
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
            let t = std::time::Instant::now();
            self.query_stats.borrow_mut().t_invert_pos +=
                t.duration_since(last).subsec_nanos() as usize;
            last = t;

            let start = plain_pos.checked_sub(offset).or_else(|| {
                self.query_stats.borrow_mut().out_of_bounds += 1;
                None
            })?;

            // The `get` fails when the minimizer match is too close to the end
            // and the pattern doesn't fit after it.
            let end = start + pattern.len();
            if end > self.seq.len() {
                self.query_stats.borrow_mut().out_of_bounds += 1;
                return None;
            }

            let matches = self.seq.slice(start..end) == pattern;
            let t = std::time::Instant::now();
            self.query_stats.borrow_mut().t_check += t.duration_since(last).subsec_nanos() as usize;
            last = t;
            if !matches {
                self.query_stats.borrow_mut().mismatches += 1;
                return None;
            }

            // Check that the range is fully inside a single input read.
            // self.ranges.index_of(start);
            let range_end = unsafe { self.ranges.succ_unchecked::<true>(start).1 };

            let t = std::time::Instant::now();
            self.query_stats.borrow_mut().t_ranges +=
                t.duration_since(last).subsec_nanos() as usize;
            last = t;

            if end > range_end {
                self.query_stats.borrow_mut().bad_ranges += 1;
                return None;
            }

            self.query_stats.borrow_mut().matches += 1;
            Some(start)
        })))
    }
}

pub fn read_human_genome() -> AsciiSeqVec {
    *INIT_TRACE;
    let _timer = Timer::new("Reading");
    let Ok(mut reader) = needletail::parse_fastx_file("human-genome.fa") else {
        panic!("Did not find human-genome.fa. Add/symlink it to test runtime on it.");
    };
    let mut seq = AsciiSeqVec::default();
    while let Some(r) = reader.next() {
        seq.push_seq(AsciiSeq(&r.unwrap().seq()));
    }
    seq
}

#[cfg(test)]
mod test {
    use indices::DivSufSortSa;
    use mem_dbg::SizeFlags;
    use packed_seq::{AsciiSeqVec, SeqVec};
    use sketchers::{IdentityParams, MinimizerParams};
    use tracing::trace;

    use super::*;

    #[test]
    fn test_identity_simple() {
        let seq = AsciiSeq(b"ACGTACGTACGTACGT");
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa { compress: false });
        let uindex = UIndex::build(seq.to_vec(), sketcher, ms_index);
        let query = AsciiSeq(b"ACGT");
        let mut occ = uindex.query(query).unwrap().collect::<Vec<_>>();
        occ.sort();
        assert_eq!(occ, vec![0, 4, 8, 12]);
    }
    #[test]
    fn test_identity_positive() {
        let seq = AsciiSeqVec::random(1000000);
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa { compress: false });
        let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
        for _ in 0..1000 {
            let len = rand::random::<usize>() % 100;
            let pos = rand::random::<usize>() % (seq.len() - len);
            let query = seq.slice(pos..pos + len);
            let occ = uindex.query(query).unwrap().collect::<Vec<_>>();
            assert!(occ.len() > 0);
            for &pos in &occ {
                assert_eq!(seq.slice(pos..pos + len), query);
            }
        }
    }
    #[test]
    fn test_identity_negative() {
        let seq = AsciiSeqVec::random(1000000);
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa { compress: false });
        let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
        for _ in 0..1000 {
            let len = 16;
            let query = AsciiSeqVec::random(len);
            let occ = uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
            assert_eq!(occ.len(), 0);
        }
    }
    #[test]
    fn test_minspace_positive() {
        let seq = AsciiSeqVec::random(1000000);

        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa { compress: false });
        let uindex = UIndex::build(seq.clone(), sketcher, ms_index);

        for remap in [false, true] {
            for l in [1, 10, 100] {
                for k in [1, 2, 3, 4, 5, 6, 7, 8] {
                    if k > l {
                        continue;
                    }
                    let sketcher = SketcherBuilderEnum::Minimizer(MinimizerParams {
                        l,
                        k,
                        remap,
                        cacheline_ef: false,
                    });
                    let index = UIndex::build(seq.clone(), sketcher, ms_index);
                    for _ in 0..1000 {
                        let len = l + rand::random::<usize>() % 100;
                        let pos = rand::random::<usize>() % (seq.len() - len);
                        let query = seq.slice(pos..pos + len);

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
        let seq = AsciiSeqVec::random(1000000);

        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa { compress: false });
        let index = UIndex::build(seq.clone(), sketcher, ms_index);

        for remap in [false, true] {
            for l in [1, 10, 100] {
                for k in [1, 2, 3, 4, 5, 6, 7, 8] {
                    if k > l {
                        continue;
                    }
                    let sketcher = SketcherBuilderEnum::Minimizer(MinimizerParams {
                        l,
                        k,
                        remap,
                        cacheline_ef: false,
                    });
                    let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
                    for _ in 0..1000 {
                        let len = l + rand::random::<usize>() % 100;
                        let query = AsciiSeqVec::random(len);

                        let mut index_occ =
                            index.query(query.as_slice()).unwrap().collect::<Vec<_>>();
                        let mut uindex_occ =
                            uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
                        index_occ.sort();
                        uindex_occ.sort();
                        assert_eq!(index_occ, uindex_occ);
                    }
                }
            }
        }
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
                    let sketcher = SketcherBuilderEnum::Minimizer(MinimizerParams {
                        l,
                        k,
                        remap,
                        cacheline_ef: false,
                    });
                    let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa { compress: false });
                    let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
                    trace!(
                        "Index size:   {} MB",
                        uindex.mem_size(SizeFlags::default()) / (1024 * 1024)
                    );
                    timer.next("Query");
                    for _ in 0..5000 {
                        let len = l + 2 * rand::random::<usize>() % l;
                        let pos = rand::random::<usize>() % (seq.len() - len);
                        let query = seq.slice(pos..pos + len);

                        let uindex_occ = uindex.query(query).unwrap().collect::<Vec<_>>();
                        assert!(uindex_occ.contains(&pos));
                    }
                }
            }
        }
    }
}
