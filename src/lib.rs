use std::{cell::RefCell, collections::HashMap};

use indices::{IndexBuilderEnum, IndexEnum};
use mem_dbg::{MemDbg, MemSize, SizeFlags};
use packed_seq::{Seq, SeqVec};
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
// TODO: Count false positives and other events. Add tracing or return them directly.
//
// TODO: Drop the minimizer space string, and instead compare against plain-space minimizers.

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
    fn query<'s>(
        &'s self,
        pattern: &[u8],
        seq: impl Seq<'s>,
        sketcher: &impl Sketcher,
    ) -> Box<dyn Iterator<Item = usize> + 's>;
}

/// Sketch a plain sequence to minimizer space.
pub trait SketcherBuilder {
    type Sketcher: Sketcher + 'static;
    /// Take an input text, find its minimizers, and compress to the target space.
    /// Additionally log statistics to `stats`.
    fn sketch_with_stats<'s, S: Seq<'s>>(
        &self,
        seq: S,
        stats: &Stats,
    ) -> (Self::Sketcher, MsSequence);

    /// Take an input text, find its minimizers, and compress to the target space.
    fn sketch<'s>(&self, seq: impl Seq<'s>) -> (Self::Sketcher, MsSequence) {
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

    fn k(&self) -> usize;

    fn len(&self) -> usize;

    /// Take an input text, compute its minimizers, and compress those into the
    /// target `u8` alphabet. This could be done a few ways, e.g.:
    /// - concatenating the KmerVals,
    /// - using a hash function to map the KmerVals to a smaller range.
    /// Also returns the position in `seq` of the first minimizer.
    /// Returns `None` when `seq` is too short to contain a minimizer.
    fn sketch<'s>(&self, seq: impl Seq<'s>) -> Result<(MsSequence, usize), SketchError>;

    /// Take a *byte* position of a character in the minimizer space, and return its start position in the original sequence.
    /// Returns `None` when the position in the minimizer space text is not aligned with the size of the encoded minimizers.
    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize>;

    /// Return the value of the minimizer at the given position in the sketched sequence.
    fn get_ms_minimizer(&self, ms_seq: &[u8], ms_pos: usize) -> Option<usize>;

    /// Return the value of the minimizer at the given position in the sketched sequence.
    fn get_ms_minimizer_via_plaintext<'s>(&self, seq: impl Seq<'s>, ms_pos: usize)
        -> Option<usize>;
}

#[derive(MemSize, MemDbg)]
pub struct UIndex<SV: SeqVec> {
    seq: SV,
    sketcher: SketcherEnum,
    ms_index: IndexEnum,
    query_stats: RefCell<QueryStats>,
    stats: Stats,
    ranges: sux::dict::elias_fano::EfDict,
}

#[derive(MemSize, MemDbg, Default, Debug)]
struct QueryStats {
    /// The total number of queries.
    queries: usize,
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
    /// Bad ranges.
    bad_ranges: usize,
    /// Matches.
    matches: usize,

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

impl<SV: SeqVec> Drop for UIndex<SV> {
    fn drop(&mut self) {
        let QueryStats {
            mut queries,
            too_short,
            unknown_minimizer,
            misaligned_ms_pos,
            out_of_bounds,
            mismatches,
            bad_ranges,
            matches,
            mut t_sketch,
            mut t_search,
            mut t_invert_pos,
            mut t_check,
            mut t_ranges,
        } = self.query_stats.take();

        queries = queries.max(1);

        t_sketch /= queries;
        t_search /= queries;
        t_invert_pos /= queries;
        t_check /= queries;
        t_ranges /= queries;

        tracing::info!(
            "QUERY STATS:
queries           {queries:>9}
too short         {too_short:>9}
unknown minimizer {unknown_minimizer:>9}
misaligned ms pos {misaligned_ms_pos:>9}
out of bounds     {out_of_bounds:>9}
mismatches        {mismatches:>9}
bad_ranges        {bad_ranges:>9}
matches           {matches:>9}
t_sketch          {t_sketch:>9} ns/query
t_search          {t_search:>9} ns/query
t_invert_pos      {t_invert_pos:>9} ns/query
t_check           {t_check:>9} ns/query
t_ranges          {t_ranges:>9} ns/query"
        );
    }
}

impl<SV: SeqVec> UIndex<SV> {
    /// 1. Sketch input to minimizer space.
    /// 2. Build minimizer space index.
    pub fn build(
        mut seq: SV,
        sketch_params: SketcherBuilderEnum,
        index_params: IndexBuilderEnum,
    ) -> Self {
        *INIT_TRACE;
        let stats = Stats::default();
        let mut timer = Timer::new_stats("Sketch", &stats);
        let (sketcher, ms_seq) = sketch_params.sketch_with_stats(seq.as_slice(), &stats);
        timer.next("Build");
        let ms_index = index_params.build_with_stats(ms_seq.0, Sketcher::width(&sketcher), &stats);
        drop(timer);

        // Build seq ranges.
        let mut ef_ranges = sux::dict::elias_fano::EliasFanoBuilder::new(
            2 * seq.ranges().len(),
            seq.ranges().last().unwrap().1,
        );
        for (start, end) in std::mem::take(seq.ranges()) {
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
        let seq_size = uindex.seq.mem_size(SizeFlags::default()) as f32 / 1000000.;
        uindex.stats.add("seq_size_MB", seq_size);
        trace!("seq    size:   {seq_size:>8.3} MB",);
        let sketch_size = uindex.sketcher.mem_size(SizeFlags::default()) as f32 / 1000000.;
        uindex.stats.add("sketch_size_MB", sketch_size);
        trace!("Sketch size:   {sketch_size:>8.3} MB",);
        let ms_seq_size = uindex.ms_index.inner().seq_size() as f32 / 1000000.;
        uindex.stats.add("ms_seq_size_MB", ms_seq_size);
        trace!("ms-seq size:   {ms_seq_size:>8.3} MB",);
        let sa_size = uindex.ms_index.inner().sa_size() as f32 / 1000000.;
        uindex.stats.add("sa_size_MB", sa_size);
        trace!("SA     size:   {sa_size:>8.3} MB",);
        let ranges_size = uindex.ranges.mem_size(SizeFlags::default()) as f32 / 1000000.;
        uindex.stats.add("ranges_size_MB", ranges_size);
        trace!("Ranges size:   {ranges_size:>8.3} MB",);
        let total_size = sketch_size + ms_seq_size + sa_size + ranges_size;
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
        stats.set("query_bad_ranges", qs.bad_ranges);
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
        pattern: <SV as SeqVec>::Seq<'p>,
    ) -> Option<Box<dyn Iterator<Item = usize> + 'p>> {
        self.query_stats.borrow_mut().queries += 1;
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
        let ms_occ = self
            .ms_index
            .query(&ms_pattern.0, self.seq.as_slice(), &self.sketcher);
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

pub fn read_chromosomes<SV: SeqVec>(cnt_max: usize) -> SV {
    *INIT_TRACE;
    let _timer = Timer::new("Reading");
    let Ok(mut reader) = needletail::parse_fastx_file("human-genome.fa") else {
        panic!("Did not find human-genome.fa. Add/symlink it to test runtime on it.");
    };
    let mut seq = SV::default();
    let mut cnt = 0;
    while let Some(r) = reader.next() {
        seq.push_ascii(&r.unwrap().seq());
        cnt += 1;
        if cnt == cnt_max {
            break;
        }
    }
    trace!(
        "Read human genome: {cnt} chromosomes of total length {}Mbp and size {}MB",
        seq.len() / 1000000,
        seq.mem_size(SizeFlags::default()) / 1000000
    );
    seq
}

pub fn read_human_genome<SV: SeqVec>() -> SV {
    read_chromosomes(usize::MAX)
}

#[cfg(test)]
mod test {
    use indices::DivSufSortSa;
    use packed_seq::{AsciiSeqVec, PackedSeqVec, SeqVec};
    use sketchers::{IdentityParams, MinimizerParams};
    use tracing::trace;

    use super::*;

    #[test]
    fn test_identity_simple() {
        let seq = AsciiSeqVec::from_ascii(b"ACGTACGTACGTACGT");
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: true,
            compress: true,
        });
        let uindex = UIndex::build(seq, sketcher, ms_index);
        let query = AsciiSeqVec::from_ascii(b"ACGT");
        let mut occ = uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
        occ.sort();
        assert_eq!(occ, vec![0, 4, 8, 12]);
    }
    #[test]
    fn test_identity_positive() {
        let seq = AsciiSeqVec::random(1000000);
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: true,
            compress: true,
        });
        let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
        for _ in 0..100 {
            let len = rand::random::<usize>() % 100;
            let pos = rand::random::<usize>() % (seq.len() - len);
            let query = seq.slice(pos..pos + len);
            let occ = uindex.query(query).unwrap().collect::<Vec<_>>();
            assert!(occ.contains(&pos));
            for &pos in &occ {
                assert_eq!(seq.slice(pos..pos + len), query);
            }
        }
    }
    #[test]
    fn test_identity_negative() {
        let seq = AsciiSeqVec::random(1000000);
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: true,
            compress: true,
        });
        let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
        for _ in 0..100 {
            let len = 32;
            let query = AsciiSeqVec::random(len);
            let occ = uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
            assert_eq!(occ.len(), 0);
        }
    }
    #[test]
    fn test_minspace_positive() {
        let seq = PackedSeqVec::random(1000000);

        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: true,
            compress: true,
        });

        for remap in [false, true] {
            for l in [10, 100] {
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
                    for _ in 0..100 {
                        let len = l + rand::random::<usize>() % 100;
                        let pos = rand::random::<usize>() % (seq.len() - len);
                        let query = seq.slice(pos..pos + len);

                        let uindex_occ = uindex.query(query).unwrap().collect::<Vec<_>>();
                        assert!(
                            uindex_occ.contains(&pos),
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
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: true,
            compress: true,
        });
        let index = UIndex::build(seq.clone(), sketcher, ms_index);

        for remap in [false, true] {
            for l in [10, 100] {
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
                    for _ in 0..100 {
                        let len = l + rand::random::<usize>() % 100;
                        let query = AsciiSeqVec::random(len);

                        let mut index_occ =
                            index.query(query.as_slice()).unwrap().collect::<Vec<_>>();
                        let mut uindex_occ =
                            uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
                        index_occ.sort();
                        uindex_occ.sort();
                        assert_eq!(index_occ, uindex_occ, "l {l} k {k} remap {remap}");
                    }
                }
            }
        }
    }
    #[test]
    fn test_identity_positive_noms() {
        let seq = AsciiSeqVec::random(1000000);
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: false,
            compress: true,
        });
        let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
        for _ in 0..100 {
            let len = rand::random::<usize>() % 100;
            let pos = rand::random::<usize>() % (seq.len() - len);
            let query = seq.slice(pos..pos + len);
            let occ = uindex.query(query).unwrap().collect::<Vec<_>>();
            assert!(occ.contains(&pos));
            for &pos in &occ {
                assert_eq!(seq.slice(pos..pos + len), query);
            }
        }
    }
    #[test]
    fn test_identity_negative_noms() {
        let seq = AsciiSeqVec::random(1000000);
        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: false,
            compress: true,
        });
        let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
        for _ in 0..100 {
            let len = 16;
            let query = AsciiSeqVec::random(len);
            let occ = uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
            assert_eq!(occ.len(), 0);
        }
    }
    #[test]
    fn test_minspace_positive_noms() {
        let seq = PackedSeqVec::random(1000000);

        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: false,
            compress: true,
        });

        for remap in [false] {
            for l in [10, 100] {
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
                    for _ in 0..100 {
                        let len = l + rand::random::<usize>() % 100;
                        let pos = rand::random::<usize>() % (seq.len() - len);
                        let query = seq.slice(pos..pos + len);

                        let uindex_occ = uindex.query(query).unwrap().collect::<Vec<_>>();
                        assert!(
                            uindex_occ.contains(&pos),
                            "l {l} k {k} remap {remap} pos {pos} query {query:?}"
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn test_minspace_negative_noms() {
        let seq = AsciiSeqVec::random(1000000);

        let sketcher = SketcherBuilderEnum::IdentityParams(IdentityParams);
        let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: false,
            compress: true,
        });
        let index = UIndex::build(seq.clone(), sketcher, ms_index);

        for remap in [false] {
            for l in [10, 100] {
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
                    for _ in 0..100 {
                        let len = l + rand::random::<usize>() % 100;
                        let query = AsciiSeqVec::random(len);

                        let mut index_occ =
                            index.query(query.as_slice()).unwrap().collect::<Vec<_>>();
                        let mut uindex_occ =
                            uindex.query(query.as_slice()).unwrap().collect::<Vec<_>>();
                        index_occ.sort();
                        uindex_occ.sort();
                        assert_eq!(index_occ, uindex_occ, "l {l} k {k} remap {remap}");
                    }
                }
            }
        }
    }

    #[test]
    #[ignore = "needs human-genome.fa"]
    fn human_genome() {
        let seq: PackedSeqVec = read_chromosomes(1);

        let ql = 256;
        let compress = true;
        for (k, l) in [(8, 32), (16, 64)] {
            for remap in [false, true] {
                for store_seq in [false, true] {
                    let mut timer = Timer::new("Build");
                    trace!("remap {remap} l {l} k {k}");
                    let sketcher = SketcherBuilderEnum::Minimizer(MinimizerParams {
                        l,
                        k,
                        remap,
                        cacheline_ef: false,
                    });
                    let ms_index = IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
                        store_ms_seq: store_seq,
                        compress,
                    });
                    let uindex = UIndex::build(seq.clone(), sketcher, ms_index);
                    timer.next("Query");
                    for _ in 0..100000 {
                        let pos = rand::random::<usize>() % (seq.len() - ql);
                        let query = seq.slice(pos..pos + ql);

                        let bad_ranges_before = uindex.query_stats.borrow().bad_ranges;
                        let uindex_occ = uindex.query(query).unwrap().collect::<Vec<_>>();
                        let bad_ranges_after = uindex.query_stats.borrow().bad_ranges;
                        assert!(uindex_occ.contains(&pos) || bad_ranges_after != bad_ranges_before);
                    }
                }
            }
        }
    }
}
