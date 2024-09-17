mod indices;
mod sketchers;

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
// - Currently we only take `u8` alphabet as input. That should be extended to
//   2-bit encoded data.
// TODO: Generalize underlying KmerVal and Pos types.
// TODO: Streaming input
// TODO: Rolling hash for O(1) checking.
// TODO: True constant time checking?
// TODO: Consider whether `Index` should own the input.

/// A packed minimizer representation.
/// Bit width of the underlying alphabet is unspecified, and should not matter:
/// really this should only be used as a unique identifier of the kmer/minimizer
/// that is only checked for equality, not inspected.
type KmerVal = u64;
/// The start position of a minimizer.
type Pos = usize;

/// A reference to a sequence over a u8 alphabet.
type Seq<'s> = &'s [u8];
/// A sequence over a u8 alphabet.
type Sequence = Vec<u8>;
/// A reference to a minimizer space sequence.
struct MsSeq<'s>(&'s [u8]);
/// A minimizer space sequence.
pub struct MsSequence(Vec<u8>);

/// A generic index to locate strings.
/// The index owns the input text.
pub trait IndexBuilder {
    type Index: Index + 'static;
    /// Build an index on the text.
    fn build(text: Sequence) -> Self::Index;
}

pub trait Index {
    /// Return all places where the pattern occurs.
    fn query<'i>(&'i self, pattern: &[u8]) -> Box<dyn Iterator<Item = usize> + 'i>;
}

/// Sketch a plain sequence to minimizer space.
pub trait SketcherBuilder {
    type Sketcher: Sketcher + 'static;
    /// Take an input text, find its minimizers, and compress to the target space.
    fn build(k: usize, w: usize, seq: &[u8]) -> (Self::Sketcher, MsSequence);
}

pub trait Sketcher {
    /// Take an input text, compute its minimizers, and compress those into the
    /// target `u8` alphabet. This could be done a few ways, e.g.:
    /// - concatenating the KmerVals,
    /// - using a hash function to map the KmerVals to a smaller range.
    /// Also returns the position in `seq` of the first minimizer.
    /// Returns `None` when `seq` is too short to contain a minimizer.
    fn to_ms(&self, seq: &[u8]) -> Option<(MsSequence, usize)>;

    /// Take a position of a character in the minimizer space, and return its start position in the original sequence.
    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> usize;
}

pub struct UIndex<'s> {
    seq: Seq<'s>,
    sketcher: Box<dyn Sketcher>,
    ms_index: Box<dyn Index>,
}

impl<'s> UIndex<'s> {
    /// 1. Sketch input to minimizer space.
    /// 2. Build minimiser space sequence.
    /// 2. (optionally) remap to smaller alphabet
    /// 3. Concat the minimizer values.
    /// 4. Reinterpret as underlying 2bit alphabet.
    /// 5. Build underlying index.
    pub fn build<C: SketcherBuilder, I: IndexBuilder>(k: usize, w: usize, seq: Seq<'s>) -> Self {
        let (sketcher, ms_seq) = C::build(k, w, seq);
        let ms_index = Box::new(I::build(ms_seq.0));
        Self {
            seq,
            sketcher: Box::new(sketcher),
            ms_index,
        }
    }

    /// 1. Sketch query to minimizer space.
    /// 2. Query the minimizer space index.
    /// 3. Check all occurrences.
    /// Returns `None` if the pattern is too short and does not contain a minimizer.
    pub fn query<'p>(&'p self, pattern: Seq<'p>) -> Option<Box<dyn Iterator<Item = usize> + 'p>> {
        let (ms_pattern, offset) = self.sketcher.to_ms(pattern)?;
        let ms_occ = self.ms_index.query(&ms_pattern.0);
        Some(Box::new(ms_occ.filter_map(move |ms_pos| {
            // Checking:
            // 1. Map minimizer space pos back to original space.
            // 2. Check query against original sequence.
            let pos = self.sketcher.ms_pos_to_plain_pos(ms_pos) - offset;
            if &self.seq[pos..pos + pattern.len()] == pattern {
                Some(pos)
            } else {
                None
            }
        })))
    }
}
