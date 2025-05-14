//! TODO: FM-index:
//! - faster-minuter
//! - quad-wavelet-tree
#[cfg(feature = "awry")]
mod fm_awry;
#[cfg(feature = "bio")]
mod fm_bio;
mod fm_sdsl;
mod sa_divsufsort;
mod sa_libsais;
mod suffix_array;

#[cfg(feature = "awry")]
pub use fm_awry::FmAwryParams;
#[cfg(feature = "bio")]
pub use fm_bio::FmBioParams;
pub use fm_sdsl::FmSdslParams;
pub use sa_divsufsort::DivSufSortSa;
pub use sa_libsais::LibSaisSa;
