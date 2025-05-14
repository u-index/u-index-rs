//! TODO: FM-index:
//! - faster-minuter
//! - quad-wavelet-tree
#[cfg(feature = "awry")]
mod fm_awry;
#[cfg(feature = "bio")]
mod fm_bio;
#[cfg(feature = "sdsl-lite-fm")]
mod fm_sdsl;
#[cfg(feature = "libdivsufsort-rs")]
mod sa_divsufsort;
mod sa_libsais;
mod suffix_array;

#[cfg(feature = "awry")]
pub use fm_awry::FmAwryParams;
#[cfg(feature = "bio")]
pub use fm_bio::FmBioParams;
#[cfg(feature = "sdsl-lite-fm")]
pub use fm_sdsl::FmSdslParams;
#[cfg(feature = "libdivsufsort-rs")]
pub use sa_divsufsort::DivSufSortSa;
pub use sa_libsais::LibSaisSa;
