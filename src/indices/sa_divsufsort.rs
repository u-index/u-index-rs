use mem_dbg::{MemDbg, MemSize, SizeFlags};
use tracing::trace;

use crate::{
    utils::{Stats, Timer},
    Index, IndexBuilder, Sketcher, S,
};

/// Build a 32-bit suffix array using `libdivsufsort`.
#[derive(Clone, Copy)]
pub struct DivSufSortSa {
    pub store_ms_seq: bool,
    pub compress: bool,
}

impl IndexBuilder for DivSufSortSa {
    type Index = SuffixArray;

    fn build_with_stats(&self, ms_seq: Vec<u8>, width: usize, stats: &Stats) -> Self::Index {
        let mut timer = Timer::new_stats("Building suffix array", stats);
        trace!("MS sequence length {}", ms_seq.len());
        stats.set("sequence length", ms_seq.len());

        // If we do not store the ms_seq, first invert the byte order of each minimizer to make sorting aligned with sorting packed u64's.

        let mut sa = libdivsufsort_rs::divsufsort(&ms_seq).expect("suffix array");

        if self.compress {
            timer.next("Compress SA");
            sa.retain(|x| *x % width as i32 == 0);
        }

        SuffixArray {
            sa,
            ms_seq: self.store_ms_seq.then(|| ms_seq),
        }
    }
}

/// A 32-bit suffix array that owns the corresponding text.
/// Uses `libdivsufsort` for searching.
#[derive(MemSize, MemDbg)]
pub struct SuffixArray {
    ms_seq: Option<Vec<u8>>,
    sa: Vec<i32>,
}

impl SuffixArray {
    pub fn seq_size(&self) -> usize {
        self.ms_seq.mem_size(SizeFlags::default())
    }
    pub fn sa_size(&self) -> usize {
        self.sa.mem_size(SizeFlags::default())
    }

    #[inline(always)]
    fn compare_minimizers<'i>(
        &self,
        seq: S<'i>,
        p: &[u8],
        i: usize,
        j: usize,
        sketcher: &impl Sketcher,
    ) -> i32 {
        let cmp = if let Some(ms_seq) = &self.ms_seq {
            let w = sketcher.width();
            let t = &ms_seq[w * i..w * (i + 1)];
            let p = &p[w * j..w * (j + 1)];
            t.cmp(p)
        } else {
            let w = sketcher.width();
            let t = sketcher.get_ms_minimizer_via_plaintext(seq, w * i).unwrap();
            let p = sketcher.get_ms_minimizer(p, w * j).unwrap();
            t.cmp(&p)
        };
        match cmp {
            std::cmp::Ordering::Less => -1,
            std::cmp::Ordering::Equal => 0,
            std::cmp::Ordering::Greater => 1,
        }
    }

    fn compare<'i>(
        &self,
        sketcher: &impl Sketcher,
        seq: S<'i>,
        p: &[u8],
        suf: i32,
        match_: &mut i32,
    ) -> i32 {
        let w = sketcher.width() as i32;

        debug_assert_eq!(suf % w, 0);
        let mut i = suf / w + *match_;
        let mut j = *match_;
        let mut r = 0;

        let ms_seq_len = sketcher.len() as i32;
        let pattern_len = (p.len() as i32) / w;
        debug_assert!(p.len() as i32 % w == 0);
        if let Some(ms_seq) = &self.ms_seq {
            debug_assert!(ms_seq.len() as i32 % w == 0);
        }

        while i < ms_seq_len as i32 && j < pattern_len {
            r = self.compare_minimizers(seq, p, i as usize, j as usize, sketcher);
            // r = t[i as usize] as i32 - p[j as usize] as i32;
            if r != 0 {
                break;
            }
            i += 1;
            j += 1;
        }
        *match_ = j;
        if r == 0 {
            if w * j != p.len() as i32 {
                r = -1;
            }
        }
        r
    }

    // Transcribed from
    // https://github.com/y-256/libdivsufsort/blob/5f60d6f026c30fb4ac296f696b3c8b0eb71bd428/lib/utils.c
    /// Search text `t` for pattern `p` given (sparse) suffix array `sa`.
    /// Returns a `(pos, cnt)` pair where `pos` is the index of the first match and `cnt` is the number of matches.
    fn sa_search<'i>(&self, sketcher: &impl Sketcher, seq: S<'i>, p: &[u8]) -> (i32, i32) {
        let mut size = self.sa.len() as i32;
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

        if sketcher.len() == 0 || self.sa.is_empty() {
            return (0, 0);
        }
        if p.is_empty() {
            return (0, self.sa.len() as i32);
        }

        while size > 0 {
            half = size / 2;
            match_ = lmatch.min(rmatch);
            r = self.compare(
                sketcher,
                seq,
                p,
                self.sa[i as usize + half as usize],
                &mut match_,
            );
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
                    r = self.compare(
                        sketcher,
                        seq,
                        p,
                        self.sa[j as usize + half as usize],
                        &mut lmatch,
                    );
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
                    r = self.compare(
                        sketcher,
                        seq,
                        p,
                        self.sa[k as usize + half as usize],
                        &mut rmatch,
                    );
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
}

impl Index for SuffixArray {
    fn query<'i>(
        &'i self,
        pattern: &[u8],
        seq: S<'i>,
        sketcher: &impl Sketcher,
    ) -> Box<dyn Iterator<Item = usize> + 'i> {
        let (pos, cnt) = self.sa_search(sketcher, seq, pattern);
        return Box::new((pos..pos + cnt).map(move |i| self.sa[i as usize] as usize));
    }
}
