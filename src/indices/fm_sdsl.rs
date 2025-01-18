use std::{any::type_name, marker::PhantomData};

use crate::{Index, IndexBuilder};
use mem_dbg::MemSize;
use packed_seq::SeqVec;
use sdsl_lite_fm::*;
use serde_json::Value;
use tracing::{info, trace, warn};

pub struct FmSdslParams<T: SdslFmIndex<C>, C> {
    _c: PhantomData<C>,
    _t: PhantomData<T>,
}

impl<T: SdslFmIndex<C>, C> std::fmt::Debug for FmSdslParams<T, C> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct(type_name::<T>()).finish()
    }
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

impl<T: SdslFmIndex<C> + 'static, C: 'static + Copy, SV: SeqVec> IndexBuilder<SV>
    for FmSdslParams<T, C>
where
    FmSdsl<T, C>: Index<SV>,
{
    fn try_build_with_stats(
        &self,
        mut text: Vec<u8>,
        width: usize,
        stats: &crate::utils::Stats,
    ) -> Option<Box<dyn Index<SV>>> {
        stats.set_val("index", Value::String("FM-sdsl".to_string()));
        stats.set("index_width", width);
        info!("Building INT SDSL on length {}", text.len());

        // Convert from big endian to little endian.
        // SDSL has a comment saying that they want big-endian input, but their
        // implementation reads little-endian input. See
        // https://github.com/simongog/sdsl-lite/issues/418

        // If C==u8 and width>1, explicitly check for zero bytes.
        if std::mem::size_of::<C>() == 1 && width > 1 {
            if text.iter().find(|x| **x == 0).is_some() {
                warn!("Input contains a zero byte!");
                return None;
            }
        }

        // If C==u8 and width>1, explicitly check for zero bytes.
        if std::mem::size_of::<C>() == 8 && width > 4 {
            warn!("SDSL int-FM index does not like large integers.");
            return None;
        }

        // Check for zero values.
        let mut has_zero = false;
        let zero = &[0; 8][..width];
        text.chunks_mut(width).for_each(|chunk| {
            chunk.reverse();
            if chunk == zero {
                has_zero = true;
            }
        });
        if has_zero {
            return None;
        }

        let path = "/tmp/input";
        std::fs::write(path, text).unwrap();
        trace!("Written to /tmp/input");
        trace!("width: {}", width);

        Some(Box::new(FmSdsl::<T, C> {
            fm: T::new(path, width),
            _phantom_c: PhantomData,
        }))
    }
}

impl<T: SdslFmIndex<u64>, SV: SeqVec + 'static> Index<SV> for FmSdsl<T, u64> {
    fn query(
        &self,
        pattern: &[u8],
        _seq: SV::Seq<'_>,
        sketcher: &dyn crate::Sketcher<SV>,
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

impl<T: SdslFmIndex<u8>, SV: SeqVec> Index<SV> for FmSdsl<T, u8> {
    fn query(
        &self,
        pattern: &[u8],
        _seq: SV::Seq<'_>,
        _sketcher: &dyn crate::Sketcher<SV>,
    ) -> Box<dyn Iterator<Item = usize> + '_> {
        let positions = self.fm.locate(&pattern);
        let len = positions.len();
        Box::new((0..len).map(move |i| positions.get(i)))
    }
}
