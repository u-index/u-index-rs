use mem_dbg::MemSize;
use packed_seq::*;

use crate::{utils::Stats, MsSequence};

/// A generic index to locate strings.
/// The index owns the input text.
pub trait IndexBuilder<SV: SeqVec>: std::fmt::Debug {
    /// Build an index on the text, and keep track of statistics.
    /// `width` gives the width of each minimizer in bytes.
    /// Effectively, only a suffix array on the `width` byte wide integers is needed.
    fn build_with_stats(&self, text: Vec<u8>, width: usize, stats: &Stats) -> Box<dyn Index<SV>> {
        self.try_build_with_stats(text, width, stats).unwrap()
    }

    fn try_build_with_stats(
        &self,
        text: Vec<u8>,
        width: usize,
        stats: &Stats,
    ) -> Option<Box<dyn Index<SV>>> {
        Some(self.build_with_stats(text, width, stats))
    }

    /// Build an index on the text.
    fn build(&self, text: Vec<u8>, width: usize) -> Box<dyn Index<SV>> {
        self.build_with_stats(text, width, &Stats::default())
    }
}

// FIXME: Re-add MemDbg super trait.
pub trait Index<SV: SeqVec>: MemSize {
    /// Return all places where the pattern occurs.
    fn query(
        &self,
        pattern: &[u8],
        seq: SV::Seq<'_>,
        sketcher: &dyn Sketcher<SV>,
    ) -> Box<dyn Iterator<Item = usize> + '_>;
}

/// Sketch a plain sequence to minimizer space.
pub trait SketcherBuilder<SV: SeqVec>: std::fmt::Debug {
    /// Take an input text, find its minimizers, and compress to the target space.
    /// Additionally log statistics to `stats`.
    fn sketch_with_stats<'x>(
        &self,
        seq: SV::Seq<'_>,
        stats: &Stats,
    ) -> (Box<dyn Sketcher<SV>>, MsSequence);

    /// Take an input text, find its minimizers, and compress to the target space.
    fn sketch(&self, seq: SV::Seq<'_>) -> (Box<dyn Sketcher<SV>>, MsSequence) {
        self.sketch_with_stats(seq, &Stats::default())
    }
}

pub enum SketchError {
    /// The sequence is too short to be sketched, and can not be searched.
    TooShort,
    /// The sequence was sketched, but contains unknown minimizers, and thus is
    /// surely not contained in the text.
    UnknownMinimizer,
}

pub trait Sketcher<SV: SeqVec>: MemSize {
    /// Returns the width in bytes of each minimizer.
    fn width(&self) -> usize;

    fn k(&self) -> usize;

    fn len(&self) -> usize;

    /// Take an input text, compute its minimizers, and compress those into the
    /// target `u8` alphabet. This could be done a few ways, e.g.:
    /// - concatenating the KmerVals,
    /// - using a hash function to map the KmerVals to a smaller range.
    /// Also returns the position in `seq` of the first minimizer.
    /// Returns `None` when `seq` is too short to contain a minimizer.
    fn sketch(&self, seq: SV::Seq<'_>) -> Result<(MsSequence, usize), SketchError>;

    /// Take a *byte* position of a character in the minimizer space, and return its start position in the original sequence.
    /// Returns `None` when the position in the minimizer space text is not aligned with the size of the encoded minimizers.
    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize>;

    /// Return the value of the minimizer at the given position in the sketched sequence.
    fn get_ms_minimizer(&self, ms_seq: &[u8], ms_pos: usize) -> Option<usize>;

    /// Return the value of the minimizer at the given position in the sketched sequence.
    fn get_ms_minimizer_via_plaintext(&self, seq: SV::Seq<'_>, ms_pos: usize) -> Option<usize>;
}
