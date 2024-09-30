//! TODO: FM-index:
//! - faster-minuter
//! - quad-wavelet-tree
mod sa_divsufsort;

use mem_dbg::{MemDbg, MemSize};
use packed_seq::Seq;
pub use sa_divsufsort::DivSufSortSa;
use sa_divsufsort::SuffixArray;

use crate::{utils::Stats, Index, IndexBuilder, Sketcher};

#[derive(Clone, Copy)]
pub enum IndexBuilderEnum {
    DivSufSortSa(DivSufSortSa),
}

#[derive(MemSize, MemDbg)]
pub enum IndexEnum {
    SuffixArray(SuffixArray),
}

impl IndexEnum {
    pub fn inner(&self) -> &SuffixArray {
        match self {
            IndexEnum::SuffixArray(index) => &index,
        }
    }
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
    fn query<'s>(
        &'s self,
        pattern: &[u8],
        seq: impl Seq<'s>,
        sketcher: &impl Sketcher,
    ) -> Box<dyn Iterator<Item = usize> + 's> {
        match self {
            IndexEnum::SuffixArray(index) => index.query(pattern, seq, sketcher),
        }
    }
}
