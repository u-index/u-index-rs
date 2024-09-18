use mem_dbg::{MemDbg, MemSize};
use tracing::trace;

use crate::{utils::Timer, Index, IndexBuilder, Sequence};

/// Build a 32-bit suffix array using `libdivsufsort`.
#[derive(Clone, Copy)]
pub struct DivSufSortSa;

impl IndexBuilder for DivSufSortSa {
    type Index = SuffixArray;

    fn build(&self, seq: Sequence) -> Self::Index {
        let timer = Timer::new("Building suffix array");
        trace!("MS sequence length {}", seq.len());
        let sa = libdivsufsort_rs::divsufsort(&seq).expect("suffix array");
        drop(timer);

        SuffixArray { sa, seq }
    }
}

/// A 32-bit suffix array that owns the corresponding text.
/// Uses `libdivsufsort` for searching.
#[derive(MemSize, MemDbg)]
pub struct SuffixArray {
    seq: Sequence,
    sa: Vec<i32>,
}

impl Index for SuffixArray {
    fn query<'i>(&'i self, pattern: &[u8]) -> Box<dyn Iterator<Item = usize> + 'i> {
        let (pos, cnt) =
            libdivsufsort_rs::sa_search(&self.seq, pattern, &self.sa).expect("sa search");
        Box::new((pos..pos + cnt).map(move |i| self.sa[i as usize] as usize))
    }
}
