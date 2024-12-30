//! TODO: FM-index:
//! - faster-minuter
//! - quad-wavelet-tree
mod fm_awry;
mod fm_bio;
mod sa_divsufsort;
mod sa_libsais;
mod suffix_array;

pub use fm_awry::FmAwryParams;
pub use fm_bio::FmBioParams;
use mem_dbg::MemSize;
use packed_seq::Seq;
pub use sa_divsufsort::DivSufSortSa;
pub use sa_libsais::LibSaisSa;
use suffix_array::SuffixArray;

use crate::{utils::Stats, Index, IndexBuilder, Sketcher};

#[derive(Clone, Copy)]
pub enum IndexBuilderEnum {
    DivSufSortSa(DivSufSortSa),
    LibSaisSa(LibSaisSa),
    FmBio(fm_bio::FmBioParams),
    FmAwry(fm_awry::FmAwryParams),
}

#[derive(MemSize)]
pub enum IndexEnum {
    SuffixArray(SuffixArray),
    FmBio(fm_bio::FmBio),
    FmAwry(fm_awry::FmAwry),
}

impl IndexEnum {
    pub fn log_sizes(&self, stats: &Stats) {
        match self {
            IndexEnum::SuffixArray(index) => index.log_sizes(stats),
            IndexEnum::FmBio(index) => index.log_sizes(stats),
            IndexEnum::FmAwry(_index) => {}
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
            IndexBuilderEnum::LibSaisSa(builder) => {
                IndexEnum::SuffixArray(builder.build_with_stats(text, width, stats))
            }
            IndexBuilderEnum::FmBio(builder) => {
                IndexEnum::FmBio(builder.build_with_stats(text, width, stats))
            }
            IndexBuilderEnum::FmAwry(builder) => {
                IndexEnum::FmAwry(builder.build_with_stats(text, width, stats))
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
            IndexEnum::FmBio(index) => index.query(pattern, seq, sketcher),
            IndexEnum::FmAwry(index) => index.query(pattern, seq, sketcher),
        }
    }
}
