//! TODO: FM-index:
//! - faster-minuter
//! - quad-wavelet-tree
mod sa_divsufsort;

use mem_dbg::{MemDbg, MemSize};
pub use sa_divsufsort::DivSufSortSa;
use sa_divsufsort::SuffixArray;

use crate::{utils::Stats, Index, IndexBuilder};

#[derive(Clone, Copy)]
pub enum IndexBuilderEnum {
    DivSufSortSa(DivSufSortSa),
}

#[derive(MemSize, MemDbg)]
pub enum IndexEnum {
    SuffixArray(SuffixArray),
}

impl IndexBuilder for IndexBuilderEnum {
    type Index = IndexEnum;

    fn build_with_stats(&self, text: Vec<u8>, width: usize, stats: &Stats) -> Self::Index {
        match self {
            IndexBuilderEnum::DivSufSortSa(builder) => {
                IndexEnum::SuffixArray(builder.build_with_stats(text, width, stats))
            }
        }
    }
}

impl Index for IndexEnum {
    fn query<'i>(&'i self, pattern: &[u8]) -> Box<dyn Iterator<Item = usize> + 'i> {
        match self {
            IndexEnum::SuffixArray(index) => index.query(pattern),
        }
    }
}
