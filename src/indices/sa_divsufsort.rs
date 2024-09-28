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
        let (pos, cnt) =
            libdivsufsort_rs::sa_search(&self.seq, pattern, &self.sa).expect("sa search");
        Box::new((pos..pos + cnt).map(move |i| self.sa[i as usize] as usize))
    }
}
