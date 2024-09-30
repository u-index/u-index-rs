//! Exported functions:
//! 1. Read human genome.
//! 2. Build u-index either plain or in mini space.
//! 3. bench_positive queries.
//! 4. Get stats.

use super::*;
use indices::DivSufSortSa;
use pyo3::prelude::*;
use sketchers::{IdentityParams, MinimizerParams};

#[pyfunction]
pub fn read() -> PyResult<SV> {
    Ok(super::read_human_genome())
}

#[pyfunction]
pub fn build_plain(seq: SV) -> PyResult<UIndex> {
    Ok(UIndex::build(
        seq,
        SketcherBuilderEnum::IdentityParams(IdentityParams),
        IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq: false,
            compress: false,
        }),
    ))
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
) -> PyResult<UIndex> {
    Ok(UIndex::build(
        seq,
        SketcherBuilderEnum::Minimizer(MinimizerParams {
            k,
            l,
            remap,
            cacheline_ef,
        }),
        IndexBuilderEnum::DivSufSortSa(DivSufSortSa {
            store_ms_seq,
            compress,
        }),
    ))
}

#[pyfunction]
pub fn gen_queries(uindex: &UIndex, len: usize, count: usize) -> PyResult<Vec<(usize, usize)>> {
    Ok(uindex.gen_query_positions(len, count))
}

#[pyfunction]
pub fn bench(uindex: &UIndex, queries: Vec<(usize, usize)>) -> PyResult<f64> {
    Ok(uindex.bench_positive(&queries))
}

#[pyfunction]
pub fn stats(uindex: &UIndex) -> PyResult<HashMap<&'static str, f32>> {
    Ok(uindex.stats())
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn uindex(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read, m)?)?;
    m.add_function(wrap_pyfunction!(build_plain, m)?)?;
    m.add_function(wrap_pyfunction!(build_minimized, m)?)?;
    m.add_function(wrap_pyfunction!(gen_queries, m)?)?;
    m.add_function(wrap_pyfunction!(self::bench, m)?)?;
    m.add_function(wrap_pyfunction!(stats, m)?)?;
    Ok(())
}
