use std::collections::HashMap;

use indices::IndexBuilderEnum;
use pyo3::pyclass;
use sketchers::SketcherBuilderEnum;

pub mod indices;
pub mod sketchers;
pub mod traits;
mod u_index;
pub mod utils;

pub use traits::*;
pub use u_index::*;

mod bench;
#[cfg(feature = "python-bindings")]
mod py;

#[cfg(test)]
mod test;

// Terminology and variables:
// - Index: a datastructure that returns all locations where a pattern matches.
// - k-mer and minimizer are the used interchangeably.
// - k: the k-mer length of the minimizers.
// - l: the length of guaranteed preserved matches.
// - w = l - k + 1: the number of k-mers in a preserved window.
// - seq: the input sequence.
//   NOTE: Avoid 'text', always use 'seq'.
// - MS: MinSpace
// - ms_seq: the minimizer-space compressed text.
// - plain seq: the original sequence.

// TODO: Generalize the alphabet.
// - Currently the underlying index is assumed to only support `u8` alphabet.
// - Currently we only take `u8` encoded 2-bit alphabet as input. That should be extended to
//   2-bit encoded data and full 8-bit alphabets.
// TODO: Generalize underlying KmerVal and Pos types.
// TODO: Streaming input
// TODO: Rolling hash for O(1) checking.
// TODO: True constant time checking?
// TODO: Consider whether `Index` should own the input.
// TODO: minmers?
// TODO: Randomize the minimizer order
// TODO: Perfect hash the minimizers
// TODO: Reserve ms ID for sequence ends.
// TODO: Count false positives and other events. Add tracing or return them directly.
//
// TODO: Drop the minimizer space string, and instead compare against plain-space minimizers.

/// A minimizer space sequence.
pub struct MsSequence(Vec<u8>);
