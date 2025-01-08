use crate::{Index, IndexBuilder};
use mem_dbg::MemSize;
use sdsl_lite_fm::*;
use tracing::{info, trace};

#[derive(Debug)]
pub struct FmSdslParams<T: SdslFmIndex<C>, C> {
    _phantom_c: std::marker::PhantomData<C>,
    _phantom_t: std::marker::PhantomData<T>,
}

impl<T: SdslFmIndex<C>, C> Clone for FmSdslParams<T, C> {
    fn clone(&self) -> Self {
        Self {
            _phantom_c: self._phantom_c.clone(),
            _phantom_t: self._phantom_t.clone(),
        }
    }
}

impl<T: SdslFmIndex<C>, C> Copy for FmSdslParams<T, C> {}

pub struct FmSdsl<T: SdslFmIndex<C>, C> {
    fm: T,
    _phantom_c: std::marker::PhantomData<C>,
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
    type Index = FmSdsl<T, C>;

    fn build_with_stats(
        &self,
        mut text: Vec<u8>,
        width: usize,
        _stats: &crate::utils::Stats,
    ) -> Self::Index {
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

        FmSdsl::<T, C> {
            fm: T::new(path, width),
            _phantom_c: std::marker::PhantomData,
        }
    }
}

impl<T: SdslFmIndex<u64>> Index for FmSdsl<T, u64> {
    fn query<'s>(
        &'s self,
        pattern: &[u8],
        _seq: impl packed_seq::Seq<'s>,
        sketcher: &impl crate::Sketcher,
    ) -> Box<dyn Iterator<Item = usize> + 's> {
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
    fn query<'s>(
        &'s self,
        pattern: &[u8],
        _seq: impl packed_seq::Seq<'s>,
        _sketcher: &impl crate::Sketcher,
    ) -> Box<dyn Iterator<Item = usize> + 's> {
        let positions = self.fm.locate(&pattern);
        let len = positions.len();
        Box::new((0..len).map(move |i| positions.get(i)))
    }
}
