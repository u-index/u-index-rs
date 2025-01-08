//! TODO: FM-index:
//! - faster-minuter
//! - quad-wavelet-tree
mod fm_awry;
mod fm_bio;
mod fm_sdsl;
mod sa_divsufsort;
mod sa_libsais;
mod suffix_array;

pub use fm_awry::FmAwryParams;
pub use fm_bio::FmBioParams;
pub use fm_sdsl::FmSdslParams;
pub use sa_divsufsort::DivSufSortSa;
pub use sa_libsais::LibSaisSa;
