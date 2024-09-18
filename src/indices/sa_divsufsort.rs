use mem_dbg::{MemDbg, MemSize};

use crate::{Index, IndexBuilder, Sequence};

/// Build a 64-bit suffix array using `libdivsufsort`.
#[derive(Clone, Copy)]
pub struct DivSufSortSa;

impl IndexBuilder for DivSufSortSa {
    type Index = SuffixArray;

    fn build(&self, seq: Sequence) -> Self::Index {
        SuffixArray {
            sa: libdivsufsort_rs::divsufsort64(&seq).expect("suffix array"),
            seq,
        }
    }
}

/// A 64-bit suffix array that owns the corresponding text.
/// Uses `libdivsufsort` for searching.
#[derive(MemSize, MemDbg)]
pub struct SuffixArray {
    seq: Sequence,
    sa: Vec<i64>,
}

impl Index for SuffixArray {
    fn query<'i>(&'i self, pattern: &[u8]) -> Box<dyn Iterator<Item = usize> + 'i> {
        let (pos, cnt) =
            libdivsufsort_rs::sa_search64(&self.seq, pattern, &self.sa).expect("sa search");
        Box::new((pos..pos + cnt).map(move |i| self.sa[i as usize] as usize))
    }
}
