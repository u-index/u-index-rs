//! This file is a crime; mostly copied from `u_index.rs` and `suffix_array.rs`.
//! But this uses slightly different types and modifying all the trais was pain.

use std::{cell::RefCell, collections::HashMap};

use crate::{
    indices::IndexBuilderEnum,
    sketchers::{MinimizerParams, SketcherBuilderEnum},
    traits::*,
    utils::*,
    QueryStats,
};
use crate::{utils::Stats, Sketcher};
use itertools::Itertools;
use mem_dbg::{MemDbg, MemSize, SizeFlags};
use packed_seq::Seq;
use packed_seq::*;
use std::cmp::Ordering;
use sux::traits::SuccUnchecked;
use tracing::trace;

#[derive(MemSize)]
pub struct SIndex<SV: SeqVec> {
    pub(crate) seq: SV,
    k: usize,
    l: usize,
    ssa: SparseSuffixArray,
    pub(crate) query_stats: RefCell<QueryStats>,
    stats: Stats,
    ranges: sux::dict::elias_fano::EfDict,
}

impl<SV: SeqVec> Drop for SIndex<SV> {
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

impl<SV: SeqVec> SIndex<SV> {
    /// 1. Sketch input to minimizer space.
    /// 2. Build minimizer space index.
    pub fn build(mut seq: SV, k: usize, l: usize) -> Self {
        *INIT_TRACE;
        let stats = Stats::default();
        let mut timer = Timer::new_stats("Sketch", &stats);

        let minimizer_positions =
            minimizers::simd::minimizer::minimizer_simd_it::<false>(seq.as_slice(), k, l - k + 1)
                .dedup()
                .collect::<Vec<_>>();
        timer.next("Build");
        // TODO
        let ssa = SparseSuffixArray::new(seq.as_slice(), minimizer_positions);
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

        let sindex = Self {
            seq,
            k,
            l,
            ssa,
            query_stats: RefCell::new(QueryStats::default()),
            stats,
            ranges: ef_ranges.build_with_dict(),
        };
        let seq_size = sindex.seq.mem_size(SizeFlags::default()) as f32 / 1000000.;
        sindex.stats.add("seq_size_MB", seq_size);
        trace!("seq    size:   {seq_size:>8.3} MB",);

        sindex.ssa.log_sizes(&sindex.stats);

        let index_size = sindex.ssa.mem_size(SizeFlags::default()) as f32 / 1000000.;
        sindex.stats.add("index_size_MB", index_size);
        trace!("Index size:   {index_size:>8.3} MB",);

        let ranges_size = sindex.ranges.mem_size(SizeFlags::default()) as f32 / 1000000.;
        sindex.stats.add("ranges_size_MB", ranges_size);
        trace!("Ranges size:   {ranges_size:>8.3} MB",);

        let total_size = index_size + ranges_size;
        sindex.stats.add("total_size_MB", total_size);
        trace!("Total  size:   {total_size:>8.3} MB",);
        sindex
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
        // TODO: Find position of first minimizer of pattern.

        if pattern.len() < self.l {
            self.query_stats.borrow_mut().too_short += 1;
            return None;
        }

        let offset = minimizers::simd::minimizer::minimizer_window_naive::<false>(
            pattern.slice(0..self.l),
            self.k,
        );

        let t2 = std::time::Instant::now();
        self.query_stats.borrow_mut().t_sketch += t2.duration_since(t1).subsec_nanos() as usize;
        let ms_occ = self
            .ssa
            .query(self.seq.as_slice(), pattern.slice(offset..pattern.len()));
        let t3 = std::time::Instant::now();
        self.query_stats.borrow_mut().t_search += t3.duration_since(t2).subsec_nanos() as usize;
        let mut last = t3;
        Some(Box::new(ms_occ.filter_map(move |pos| {
            // Checking:
            // 1. Check first offset characters of pattern against original sequence.

            let start = pos.checked_sub(offset).or_else(|| {
                self.query_stats.borrow_mut().out_of_bounds += 1;
                None
            })?;

            // The `get` fails when the minimizer match is too close to the end
            // and the pattern doesn't fit after it.
            let end = start + pattern.len();
            assert!(end <= self.seq.len(), "Pattern extends beyond the text");

            let matches = self.seq.slice(start..start + offset) == pattern.slice(0..offset);
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

/// A 32-bit suffix array that owns the corresponding text.
#[derive(MemSize, MemDbg)]
pub struct SparseSuffixArray {
    sa: Vec<u32>,
}

impl SparseSuffixArray {
    pub fn new<'i>(text: impl Seq<'i>, mut indices: Vec<u32>) -> Self {
        let n = text.len();
        indices.sort_unstable_by_key(|idx| text.slice(*idx as usize..n));
        Self { sa: indices }
    }

    pub fn sa_size(&self) -> usize {
        self.sa.mem_size(SizeFlags::default())
    }

    pub fn log_sizes(&self, stats: &Stats) {
        let sa_size = self.sa_size() as f32 / 1000000.;
        stats.add("sa_size_MB", sa_size);
        trace!("SA     size:   {sa_size:>8.3} MB",);
    }

    fn compare<'i, S: Seq<'i>>(
        &self,
        seq: S,
        p: S,
        // Byte-position in the sketched text that we compare against.
        // Must be a multiple of the kmer width.
        suf: usize,
        // Number of bytes of pattern already matched.
        // Must be a multiple of the kmer width when matching full minimizers at a time.
        match_: &mut usize,
    ) -> Ordering {
        let seq = seq.slice(suf + *match_..seq.len());
        let pat = p.slice(*match_..p.len());
        let (r, lcp) = seq.cmp_lcp(&pat);
        *match_ += lcp;
        r
    }

    // Transcribed from
    // https://github.com/y-256/libdivsufsort/blob/5f60d6f026c30fb4ac296f696b3c8b0eb71bd428/lib/utils.c
    /// Search text `t` for pattern `p` given (sparse) suffix array `sa`.
    /// Returns a `(pos, cnt)` pair where `pos` is the index of the first match and `cnt` is the number of matches.
    fn sa_search<'i, S: Seq<'i>>(&self, seq: S, p: S) -> (usize, usize) {
        let mut size = self.sa.len();
        let mut half;
        let mut match_;
        let mut lmatch = 0;
        let mut rmatch = 0;
        let mut llmatch;
        let mut lrmatch;
        let mut rlmatch;
        let mut rrmatch;
        let mut i = 0;
        let mut j = 0;
        let mut k = 0;
        let mut r;

        if p.len() == 0 {
            return (0, self.sa.len());
        }

        while size > 0 {
            half = size / 2;
            match_ = lmatch.min(rmatch);
            r = self.compare(
                seq,
                p,
                self.sa[i as usize + half as usize] as usize,
                &mut match_,
            );
            if r.is_lt() {
                i += half + 1;
                half -= (size & 1) ^ 1;
                lmatch = match_;
            } else if r.is_gt() {
                rmatch = match_;
            } else {
                let mut lsize = half;
                j = i;
                let mut rsize = size - half - 1;
                k = i + half + 1;

                /* left part */
                llmatch = lmatch;
                lrmatch = match_;
                while lsize > 0 {
                    half = lsize >> 1;
                    lmatch = llmatch.min(lrmatch);
                    r = self.compare(
                        seq,
                        p,
                        self.sa[j as usize + half as usize] as usize,
                        &mut lmatch,
                    );
                    if r.is_lt() {
                        j += half + 1;
                        half -= (lsize & 1) ^ 1;
                        llmatch = lmatch;
                    } else {
                        lrmatch = lmatch;
                    }
                    lsize = half;
                }

                /* right part */
                rlmatch = match_;
                rrmatch = rmatch;
                while rsize > 0 {
                    half = rsize >> 1;
                    rmatch = rlmatch.min(rrmatch);
                    r = self.compare(
                        seq,
                        p,
                        self.sa[k as usize + half as usize] as usize,
                        &mut rmatch,
                    );
                    if r.is_le() {
                        k += half + 1;
                        half -= (rsize & 1) ^ 1;
                        rlmatch = rmatch;
                    } else {
                        rrmatch = rmatch;
                    }
                    rsize = half;
                }
                break;
            }
            size = half;
        }

        let pos = if k - j > 0 { j } else { i };
        let cnt = k - j;
        (pos, cnt)
    }

    fn query<'s, S: Seq<'s>>(&'s self, seq: S, pattern: S) -> impl Iterator<Item = usize> + 's {
        let (pos, cnt) = self.sa_search(seq, pattern);
        (pos..pos + cnt).map(move |i| self.sa[i as usize] as usize)
    }
}
