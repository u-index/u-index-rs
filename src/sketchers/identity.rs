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
    fn sketch_with_stats(&self, seq: PackedSeq, stats: &Stats) -> (Box<dyn Sketcher>, MsSequence) {
        stats.set_val("sketcher", Value::String("identity".to_string()));
        stats.set("sketch_skip_zero", self.skip_zero as u64);
        let seq = seq
            .iter_bp()
            .collect::<Vec<_>>();
        (
            Box::new(Identity {
                params: *self,
                len: seq.len(),
            }),
            MsSequence(seq),
        )
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
        let seq = seq
            .iter_bp()
            .collect::<Vec<_>>();
        Ok((MsSequence(seq), 0))
    }

    fn ms_pos_to_plain_pos(&self, ms_pos: usize) -> Option<usize> {
        Some(ms_pos)
    }

    fn get_ms_minimizer(&self, ms_seq: &[u8], ms_pos: usize) -> Option<usize> {
        Some(ms_seq[ms_pos] as usize)
    }

    fn get_ms_minimizer_via_plaintext(&self, seq: PackedSeq, ms_pos: usize) -> Option<usize> {
        Some(seq.get(ms_pos) as usize)
    }
}
