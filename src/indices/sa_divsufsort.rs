use mem_dbg::{MemDbg, MemSize};
use tracing::trace;

use crate::{
    utils::{Stats, Timer},
    Index, IndexBuilder,
};

/// Build a 32-bit suffix array using `libdivsufsort`.
#[derive(Clone, Copy)]
pub struct DivSufSortSa {
    pub compress: bool,
}

impl IndexBuilder for DivSufSortSa {
    type Index = SuffixArray;

    fn build_with_stats(&self, seq: Vec<u8>, width: usize, stats: &Stats) -> Self::Index {
        let mut timer = Timer::new_stats("Building suffix array", stats);
        trace!("MS sequence length {}", seq.len());
        stats.set("sequence length", seq.len());
        let mut sa = libdivsufsort_rs::divsufsort(&seq).expect("suffix array");

        if self.compress {
            timer.next("Compress SA");
            sa.retain(|x| *x % width as i32 == 0);
        }

        SuffixArray { sa, seq }
    }
}

/// A 32-bit suffix array that owns the corresponding text.
/// Uses `libdivsufsort` for searching.
#[derive(MemSize, MemDbg)]
pub struct SuffixArray {
    seq: Vec<u8>,
    sa: Vec<i32>,
}

impl Index for SuffixArray {
    fn query<'i>(&'i self, pattern: &[u8]) -> Box<dyn Iterator<Item = usize> + 'i> {
        let (pos, cnt) = sa_search(&self.ms_str, pattern, &self.sa);
        Box::new((pos..pos + cnt).map(move |i| self.sa[i as usize] as usize))
    }
fn compare(t: &[u8], p: &[u8], suf: i32, match_: &mut i32) -> i32 {
    let mut i = suf + *match_;
    let mut j = *match_;
    let mut r = 0;
    while i < t.len() as i32 && j < p.len() as i32 {
        r = t[i as usize] as i32 - p[j as usize] as i32;
        if r != 0 {
            break;
        }
        i += 1;
        j += 1;
    }
    *match_ = j;
    if r == 0 {
        if j != p.len() as i32 {
            return -1;
        } else {
            return 0;
        }
    }
    r
}

// Transcribed from
// https://github.com/y-256/libdivsufsort/blob/5f60d6f026c30fb4ac296f696b3c8b0eb71bd428/lib/utils.c
/// Search text `t` for pattern `p` given (sparse) suffix array `sa`.
/// Returns a `(pos, cnt)` pair where `pos` is the index of the first match and `cnt` is the number of matches.
fn sa_search(t: &[u8], p: &[u8], sa: &[i32]) -> (i32, i32) {
    let mut size = sa.len() as i32;
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

    if t.is_empty() || sa.is_empty() {
        return (0, 0);
    }
    if p.is_empty() {
        return (0, sa.len() as i32);
    }

    while size > 0 {
        half = size / 2;
        match_ = lmatch.min(rmatch);
        r = compare(t, p, sa[i as usize + half as usize], &mut match_);
        if r < 0 {
            i += half + 1;
            half -= (size & 1) ^ 1;
            lmatch = match_;
        } else if r > 0 {
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
                r = compare(t, p, sa[j as usize + half as usize], &mut lmatch);
                if r < 0 {
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
                r = compare(t, p, sa[k as usize + half as usize], &mut rmatch);
                if r <= 0 {
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

    let idx = if k - j > 0 { j } else { i };
    (idx, k - j)
}
