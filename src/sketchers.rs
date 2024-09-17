use crate::{MsSequence, Sketcher, SketcherBuilder};

pub struct Identity;

impl SketcherBuilder for Identity {
    type Sketcher = Self;

    fn build(_k: usize, _w: usize, seq: &[u8]) -> (Self::Sketcher, crate::MsSequence) {
        (Self, MsSequence(seq.to_vec()))
    }
}

impl Sketcher for Identity {
    fn to_ms(&self, seq: &[u8]) -> Option<(MsSequence, usize)> {
        Some((MsSequence(seq.to_vec()), 0))
    }

    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> usize {
        ms_pos
    }
}
