use crate::{Index, IndexBuilder, Sequence};

/// Build a 64-bit suffix array using `libdivsufsort`.
#[derive(Clone, Copy)]
pub struct DivSufSortSaBuilder;

impl IndexBuilder for DivSufSortSaBuilder {
    type Index = SuffixArray;

    fn build(seq: Sequence) -> Self::Index {
        SuffixArray {
            sa: libdivsufsort_rs::divsufsort64(&seq).expect("suffix array"),
            seq,
        }
    }
}

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
