use packed_seq::{PackedSeq, SeqVec};

use super::*;

/// 'Sketch' the sequence to itself.
/// Convenient for testing purposes.
#[derive(Clone, Copy)]
pub struct IdentityParams;

#[derive(MemSize, MemDbg)]
pub struct Identity {
    len: usize,
}

impl SketcherBuilder for IdentityParams {
    fn sketch_with_stats(&self, seq: PackedSeq, _stats: &Stats) -> (Box<dyn Sketcher>, MsSequence) {
        assert!(
            PackedSeq::BASES_PER_BYTE == 1,
            "Identity sketcher does not work for packed sequences."
        );
        let seq = seq.to_vec().into_raw();
        (Box::new(Identity { len: seq.len() }), MsSequence(seq))
    }
}

impl Sketcher for Identity {
    fn width(&self) -> usize {
        1
    }

    fn k(&self) -> usize {
        1
    }

    fn len(&self) -> usize {
        self.len
    }

    fn sketch(&self, seq: PackedSeq) -> Result<(MsSequence, usize), SketchError> {
        Ok((MsSequence(seq.to_vec().into_raw()), 0))
    }

    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize> {
        Some(ms_pos)
    }

    fn get_ms_minimizer(&self, ms_seq: &[u8], ms_pos: usize) -> Option<usize> {
        Some(ms_seq[ms_pos] as usize)
    }

    fn get_ms_minimizer_via_plaintext(&self, seq: PackedSeq, ms_pos: usize) -> Option<usize> {
        Some(seq.get_ascii(ms_pos) as usize)
    }
}
