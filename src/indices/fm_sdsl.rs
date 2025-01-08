use std::marker::PhantomData;

use crate::{Index, IndexBuilder};
use mem_dbg::MemSize;
use packed_seq::PackedSeq;
use sdsl_lite_fm::*;
use tracing::{info, trace};

#[derive(Debug)]
pub struct FmSdslParams<T: SdslFmIndex<C>, C> {
    _c: PhantomData<C>,
    _t: PhantomData<T>,
}

impl<T: SdslFmIndex<C>, C> FmSdslParams<T, C> {
    pub fn new() -> Self {
        Self {
            _c: PhantomData,
            _t: PhantomData,
        }
    }
}

impl<T: SdslFmIndex<C>, C> Clone for FmSdslParams<T, C> {
    fn clone(&self) -> Self {
        Self {
            _c: self._c.clone(),
            _t: self._t.clone(),
        }
    }
}

impl<T: SdslFmIndex<C>, C> Copy for FmSdslParams<T, C> {}

pub struct FmSdsl<T: SdslFmIndex<C>, C> {
    fm: T,
    _phantom_c: PhantomData<C>,
}

impl<T: SdslFmIndex<C>, C> MemSize for FmSdsl<T, C> {
    fn mem_size(&self, _flags: mem_dbg::SizeFlags) -> usize {
        self.fm.size()
    }
}

impl<T: SdslFmIndex<C> + 'static, C: 'static + Copy> IndexBuilder for FmSdslParams<T, C>
where
    FmSdsl<T, C>: Index,
{
    fn build_with_stats(
        &self,
        mut text: Vec<u8>,
        width: usize,
        _stats: &crate::utils::Stats,
    ) -> Box<dyn Index> {
        info!("Build INT SDSL on length {}", text.len());

        // Convert from big endian to little endian.
        // SDSL has a comment saying that they want big-endian input, but their
        // implementation reads little-endian input. See
        // https://github.com/simongog/sdsl-lite/issues/418

        text.chunks_mut(width).for_each(|chunk| {
            chunk.reverse();
        });

        let path = "/tmp/input";
        std::fs::write(path, text).unwrap();
        trace!("Written to /tmp/input");
        trace!("width: {}", width);

        Box::new(FmSdsl::<T, C> {
            fm: T::new(path, width),
            _phantom_c: PhantomData,
        })
    }
}

impl<T: SdslFmIndex<u64>> Index for FmSdsl<T, u64> {
    fn query(
        &self,
        pattern: &[u8],
        _seq: PackedSeq,
        sketcher: &dyn crate::Sketcher,
    ) -> Box<dyn Iterator<Item = usize> + '_> {
        // Convert pattern to a [u64].
        let width = sketcher.width();
        assert!(width <= 8);

        let mut ints = vec![0u64; pattern.len() / width];
        assert_eq!(ints.len() * width, pattern.len());
        for i in 0..ints.len() {
            // convert slice to u64.
            let mut bytes = [0u8; 8];
            bytes[8 - width..].copy_from_slice(&pattern[i * width..i * width + width]);
            ints[i] = u64::from_be_bytes(bytes);
            assert!(ints[i] < 1 << (width * 8));
        }

        let positions = self.fm.locate(&ints);
        let len = positions.len();
        Box::new((0..len).map(move |i| positions.get(i) * width))
    }
}

impl<T: SdslFmIndex<u8>> Index for FmSdsl<T, u8> {
    fn query(
        &self,
        pattern: &[u8],
        _seq: PackedSeq,
        _sketcher: &dyn crate::Sketcher,
    ) -> Box<dyn Iterator<Item = usize> + '_> {
        let positions = self.fm.locate(&pattern);
        let len = positions.len();
        Box::new((0..len).map(move |i| positions.get(i)))
    }
}
