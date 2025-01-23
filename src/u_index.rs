use std::{cell::RefCell, collections::HashMap, ops::Range};

use mem_dbg::{MemDbg, MemSize, SizeFlags};
use packed_seq::*;
use serde_json::Value;
use sux::traits::SuccUnchecked;
use tracing::{info, trace};

use crate::{traits::*, utils::*};

#[derive(MemSize)]
pub struct UIndex<'s, SV: SeqVec> {
    pub(crate) seq: &'s SV,
    sketcher: Box<dyn Sketcher<SV>>,
    ms_index: Box<dyn Index<SV>>,
    pub(crate) query_stats: RefCell<QueryStats>,
    stats: Stats,
    ranges: sux::dict::elias_fano::EfDict,
}

#[derive(MemSize, MemDbg, Default, Debug)]
pub struct QueryStats {
    /// The total number of queries.
    pub queries: usize,
    /// Pattern is too short to sketch.
    pub too_short: usize,
    /// Pattern contained a minimizer that is not in the input text.
    pub unknown_minimizer: usize,
    /// A minimizer-space match was found that does not align to the minimizer alphabet size.
    pub misaligned_ms_pos: usize,
    /// A minimizer-space match was found, but its containing sequence is too short to contain the pattern.
    pub out_of_bounds: usize,
    /// Mismatch in sequence space.
    pub mismatches: usize,
    /// Bad ranges.
    pub bad_ranges: usize,
    /// Matches.
    pub matches: usize,

    /// Total time in ns of sketching queries.
    pub t_sketch: usize,
    /// Total time in ns of looking up queries in the minimizer space index.
    pub t_search: usize,
    /// Total time in ns of converting the ms position back to the plain position.
    pub t_invert_pos: usize,
    /// Total time in ns checking matches.
    pub t_check: usize,
    /// Total time in ns checking ranges.
    pub t_ranges: usize,
}

impl<'s, SV: SeqVec> Drop for UIndex<'s, SV> {
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

        tracing::info!("mismatches:  {mismatches:>9}");
        tracing::info!("matches:     {matches:>9}");
        tracing::trace!(
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

impl<'s, SV: SeqVec + 'static> UIndex<'s, SV> {
    pub fn build(
        seq: &'s SV,
        sketch_params: &dyn SketcherBuilder<SV>,
        index_params: &dyn IndexBuilder<SV>,
    ) -> Self {
        let ranges = [0..seq.len()];
        Self::build_with_ranges(seq, &ranges, sketch_params, index_params)
    }

    pub fn build_with_ranges(
        seq: &'s SV,
        ranges: &[Range<usize>],
        sketch_params: &dyn SketcherBuilder<SV>,
        index_params: &dyn IndexBuilder<SV>,
    ) -> Self {
        Self::try_build_with_ranges(seq, ranges, sketch_params, index_params)
            .expect("Failed to build UIndex")
    }

    /// 1. Sketch input to minimizer space.
    /// 2. Build minimizer space index.
    pub fn try_build_with_ranges(
        seq: &'s SV,
        ranges: &[Range<usize>],
        sketch_params: &dyn SketcherBuilder<SV>,
        index_params: &dyn IndexBuilder<SV>,
    ) -> Option<Self> {
        *INIT_TRACE;
        let stats = Stats::default();
        let mut timer = Timer::new_stats("Sketch", &stats).info();
        let (sketcher, ms_seq) = sketch_params.sketch_with_stats(seq.as_slice(), &stats);
        timer.next("Build");
        let ms_index =
            index_params.try_build_with_stats(ms_seq.0, Sketcher::width(&*sketcher), &stats)?;
        drop(timer);

        // Build seq ranges.
        let mut ef_ranges = sux::dict::elias_fano::EliasFanoBuilder::new(
            2 * ranges.len(),
            ranges.last().unwrap().end,
        );
        for r in ranges {
            ef_ranges.push(r.start);
            ef_ranges.push(r.end);
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
        info!("seq    size:    {seq_size:>8.3} MB",);

        let sketch_size = uindex.sketcher.mem_size(SizeFlags::default()) as f32 / 1000000.;
        uindex.stats.add("sketch_size_MB", sketch_size);
        info!("Sketch size:    {sketch_size:>8.3} MB",);

        // uindex.ms_index.log_sizes(&uindex.stats);

        let index_size = uindex.ms_index.mem_size(SizeFlags::default()) as f32 / 1000000.;
        uindex.stats.add("index_size_MB", index_size);
        info!("Index size:     {index_size:>8.3} MB",);

        let ranges_size = uindex.ranges.mem_size(SizeFlags::default()) as f32 / 1000000.;
        uindex.stats.add("ranges_size_MB", ranges_size);
        trace!("Ranges size:   {ranges_size:>8.3} MB",);

        let total_size = sketch_size + index_size + ranges_size;
        uindex.stats.add("total_size_MB", total_size);
        info!("Total  size:    {total_size:>8.3} MB",);
        Some(uindex)
    }

    pub fn stats(&self) -> HashMap<&'static str, Value> {
        let stats = self.stats.clone();
        let qs = self.query_stats.borrow();
        stats.set("queries", qs.queries);
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
        pattern: SV::Seq<'p>,
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
            .query(&ms_pattern.0, self.seq.as_slice(), &*self.sketcher);
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
