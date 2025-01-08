//! Exported functions:
//! 1. Read human genome.
//! 2. Build u-index either plain or in mini space.
//! 3. bench_positive queries.
//! 4. Get stats.

use crate::bench::gen_query_positions;

use super::*;
use indices::{DivSufSortSa, FmAwryParams, FmBioParams};
use packed_seq::{PackedSeqVec, SeqVec};
use pyo3::prelude::*;
use s_index::SIndex;
use sketchers::{IdentityParams, MinimizerParams};

/// For python bindings, we always use the packed version.
type SV = PackedSeqVec;

#[pyclass]
pub struct PyUIndex(UIndex<SV>);

#[pyclass]
pub struct PySIndex(SIndex<SV>);

#[pyfunction]
pub fn read() -> PyResult<SV> {
    Ok(utils::read_human_genome())
}

#[pyfunction]
pub fn build_plain(seq: SV) -> PyResult<PyUIndex> {
    Ok(PyUIndex(UIndex::build(
        seq,
        SketcherBuilderEnum::IdentityParams(IdentityParams),
        IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: false,
            compress: false,
        }),
    )))
}

#[pyfunction]
pub fn build_sindex(seq: SV, k: usize, l: usize) -> PyResult<PySIndex> {
    Ok(PySIndex(SIndex::build(seq, k, l)))
}

#[pyfunction]
pub fn build_minimized(
    seq: SV,
    k: usize,
    l: usize,
    remap: bool,
    store_ms_seq: bool,
    compress: bool,
    cacheline_ef: bool,
) -> PyResult<PyUIndex> {
    Ok(PyUIndex(UIndex::build(
        seq,
        SketcherBuilderEnum::Minimizer(MinimizerParams {
            k,
            l,
            remap,
            cacheline_ef,
            skip_zero: false,
        }),
        IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq,
            compress,
        }),
    )))
}

#[pyfunction]
pub fn build_plain_fm(seq: SV, sa_sampling: usize) -> PyResult<PyUIndex> {
    Ok(PyUIndex(UIndex::build(
        seq,
        SketcherBuilderEnum::IdentityParams(IdentityParams),
        IndexBuilderEnum::FmAwry(FmAwryParams { sa_sampling }),
    )))
}

#[pyfunction]
pub fn build_minimized_fm(
    seq: SV,
    k: usize,
    l: usize,
    remap: bool,
    cacheline_ef: bool,
    sa_sampling: usize,
) -> PyResult<PyUIndex> {
    Ok(PyUIndex(UIndex::build(
        seq,
        SketcherBuilderEnum::Minimizer(MinimizerParams {
            k,
            l,
            remap,
            cacheline_ef,
            skip_zero: false,
        }),
        IndexBuilderEnum::FmAwry(FmAwryParams { sa_sampling }),
    )))
}

#[pyfunction]
pub fn gen_queries(seq: SV, len: usize, count: usize) -> PyResult<Vec<(usize, usize)>> {
    Ok(gen_query_positions(seq.as_slice(), len, count))
}

#[pyfunction]
pub fn bench(uindex: &PyUIndex, queries: Vec<(usize, usize)>) -> PyResult<f64> {
    Ok(uindex.0.bench_positive(&queries))
}

#[pyfunction]
pub fn bench_sindex(sindex: &PySIndex, queries: Vec<(usize, usize)>) -> PyResult<f64> {
    Ok(sindex.0.bench_positive(&queries))
}

#[pyfunction]
pub fn stats(uindex: &PyUIndex) -> PyResult<HashMap<&'static str, f32>> {
    Ok(uindex.0.stats())
}

#[pyfunction]
pub fn stats_sindex(sindex: &PySIndex) -> PyResult<HashMap<&'static str, f32>> {
    Ok(sindex.0.stats())
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn uindex(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read, m)?)?;
    m.add_function(wrap_pyfunction!(build_plain, m)?)?;
    m.add_function(wrap_pyfunction!(build_sindex, m)?)?;
    m.add_function(wrap_pyfunction!(build_minimized, m)?)?;
    m.add_function(wrap_pyfunction!(build_plain_fm, m)?)?;
    m.add_function(wrap_pyfunction!(build_minimized_fm, m)?)?;
    m.add_function(wrap_pyfunction!(gen_queries, m)?)?;
    m.add_function(wrap_pyfunction!(self::bench, m)?)?;
    m.add_function(wrap_pyfunction!(self::bench_sindex, m)?)?;
    m.add_function(wrap_pyfunction!(stats, m)?)?;
    m.add_function(wrap_pyfunction!(stats_sindex, m)?)?;

    Ok(())
}
