mod minimizers;

pub use minimizers::MinimizerParams;

use crate::{MsSequence, Seq, Sketcher, SketcherBuilder};

/// 'Sketch' the sequence to itself.
/// Convenient for testing purposes.
pub struct Identity;

impl SketcherBuilder for Identity {
    type Sketcher = Self;

    fn sketch(&self, seq: Seq) -> (Self::Sketcher, crate::MsSequence) {
        (Self, MsSequence(seq.to_vec()))
    }
}

impl Sketcher for Identity {
    fn sketch(&self, seq: Seq) -> Option<(MsSequence, usize)> {
        Some((MsSequence(seq.to_vec()), 0))
    }

    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize> {
        Some(ms_pos)
    }
}
