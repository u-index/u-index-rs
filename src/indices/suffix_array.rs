use std::cmp::Ordering;

use mem_dbg::{MemDbg, MemSize, SizeFlags};
use packed_seq::Seq;

use crate::{Index, Sketcher};

/// A 32-bit suffix array that owns the corresponding text.
/// Uses `libdivsufsort` for searching.
#[derive(MemSize, MemDbg)]
pub struct SuffixArray {
    pub(super) ms_seq: Option<Vec<u8>>,
    pub(super) sa: Vec<i32>,
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
        seq: impl Seq<'i>,
        pattern: &[u8],
        // Byte index in ms_seq.
        i: usize,
        // Byte index in pattern.
        j: usize,
        sketcher: &impl Sketcher,
    ) -> Ordering {
        let w = sketcher.width() as usize;
        if let Some(ms_seq) = &self.ms_seq {
            let t = &ms_seq[i..i + w];
            let p = &pattern[j..j + w];
            t.cmp(p)
        } else {
            let t = sketcher.get_ms_minimizer_via_plaintext(seq, i).unwrap();
            let p = sketcher.get_ms_minimizer(pattern, j).unwrap();
            t.cmp(&p)
        }
    }

    fn compare<'i>(
        &self,
        sketcher: &impl Sketcher,
        seq: impl Seq<'i>,
        p: &[u8],
        // Byte-position in the sketched text that we compare against.
        // Must be a multiple of the kmer width.
        suf: i32,
        // Number of bytes of pattern already matched.
        // Must be a multiple of the kmer width when matching full minimizers at a time.
        match_: &mut i32,
    ) -> Ordering {
        let w = sketcher.width() as i32;

        debug_assert_eq!(suf % w, 0);
        debug_assert!(p.len() as i32 % w == 0);
        if let Some(ms_seq) = &self.ms_seq {
            debug_assert!(ms_seq.len() as i32 % w == 0);
        }

        let mut i = suf + *match_;
        let mut j = *match_;
        let mut r = Ordering::Equal;

        let ms_seq_len = sketcher.len() as i32 * w;
        let pattern_len = p.len() as i32;

        while i < ms_seq_len as i32 && j < pattern_len {
            r = self.compare_minimizers(seq, p, i as usize, j as usize, sketcher);
            if r != Ordering::Equal {
                break;
            }
            i += w;
            j += w;
        }
        *match_ = j;
        if r.is_eq() {
            if j != p.len() as i32 {
                r = Ordering::Less;
            }
        }
        r
    }

    // Transcribed from
    // https://github.com/y-256/libdivsufsort/blob/5f60d6f026c30fb4ac296f696b3c8b0eb71bd428/lib/utils.c
    /// Search text `t` for pattern `p` given (sparse) suffix array `sa`.
    /// Returns a `(pos, cnt)` pair where `pos` is the index of the first match and `cnt` is the number of matches.
    fn sa_search<'i>(&self, sketcher: &impl Sketcher, seq: impl Seq<'i>, p: &[u8]) -> (i32, i32) {
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
                        sketcher,
                        seq,
                        p,
                        self.sa[j as usize + half as usize],
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
                        sketcher,
                        seq,
                        p,
                        self.sa[k as usize + half as usize],
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

        let idx = if k - j > 0 { j } else { i };
        (idx, k - j)
    }
}

impl Index for SuffixArray {
    fn query<'s>(
        &'s self,
        pattern: &[u8],
        seq: impl Seq<'s>,
        sketcher: &impl Sketcher,
    ) -> Box<dyn Iterator<Item = usize> + 's> {
        let (pos, cnt) = self.sa_search(sketcher, seq, pattern);
        return Box::new((pos..pos + cnt).map(move |i| self.sa[i as usize] as usize));
    }
}
