use std::{any::type_name, collections::HashMap, ops::Range};

use packed_seq::{AsciiSeqVec, PackedSeqVec, SeqVec};
use sdsl_lite_fm::*;
use serde_json::{Number, Value};
use uindex::{
    bench::gen_query_positions,
    indices::{DivSufSortSa, FmAwryParams, FmSdslParams, LibSaisSa},
    s_index::SIndex,
    sketchers::{IdentityParams, MinimizerParams},
    utils::{read_chromosomes, Timer, INIT_TRACE},
    IndexBuilder, SketcherBuilder, UIndex,
};

fn main() {
    color_backtrace::install();
    *INIT_TRACE;

    let query_length = 128;
    let num_queries = 10000;
    // let (seq, _ranges) = read_chromosomes::<PackedSeqVec>(1);
    let seq = &std::fs::read("english.200MB").unwrap();
    let queries = gen_query_positions(seq.len(), query_length, num_queries);

    let mut all_stats = vec![];
    // run::<PackedSeqVec>(&mut all_stats, seq, query_length, &queries);
    // run::<AsciiSeqVec>(&mut all_stats, seq, query_length, &queries);
    run::<Vec<u8>>(&mut all_stats, seq, query_length, &queries);

    // Write all_stats.
    let stats_string = serde_json::to_string(&all_stats).unwrap();
    std::fs::write("stats.json", stats_string).unwrap();
}

fn run<SV: SeqVec>(
    all_stats: &mut Vec<HashMap<&str, Value>>,
    seq: &[u8],
    query_length: usize,
    queries: &Vec<(usize, usize)>,
) {
    tracing::warn!("{}", type_name::<SV>());
    let ranges = vec![0..seq.len()];
    let seq = SV::from_ascii(seq);

    for (k, l) in [
        (0, 0),
        (2, 16),
        (2, 20),
        (3, 20),
        (3, 24),
        (4, 28),
        (4, 32),
        (8, 56),
        (8, 60),
    ] {
        // SKETCHERS
        let id = &IdentityParams { skip_zero: false };
        let id_skip = &IdentityParams { skip_zero: true };
        let min_no_remap = &MinimizerParams {
            k,
            l,
            remap: false,
            cacheline_ef: true,
            skip_zero: false,
        };
        let min_no_remap_skip = &MinimizerParams {
            k,
            l,
            remap: false,
            cacheline_ef: true,
            skip_zero: true,
        };
        let min_remap = &MinimizerParams {
            k,
            l,
            remap: true,
            cacheline_ef: true,
            skip_zero: false,
        };
        let min_remap_skip = &MinimizerParams {
            k,
            l,
            remap: true,
            cacheline_ef: true,
            skip_zero: true,
        };

        // INDICES
        let sais_ms = &LibSaisSa {
            store_ms_seq: true,
            par: false,
        };
        let sais_no_ms = &LibSaisSa {
            store_ms_seq: false,
            par: false,
        };
        let awry32 = &FmAwryParams { sa_sampling: 32 };
        let awry64 = &FmAwryParams { sa_sampling: 64 };
        let sdsl_byte_32 = &FmSdslParams::<FmIndexByte32Ptr, _>::new();
        let sdsl_byte_64 = &FmSdslParams::<FmIndexByte64Ptr, _>::new();
        let sdsl_int_32 = &FmSdslParams::<FmIndexInt32Ptr, _>::new();
        let sdsl_int_64 = &FmSdslParams::<FmIndexInt64Ptr, _>::new();

        let params: Vec<(&dyn IndexBuilder<SV>, &dyn SketcherBuilder<SV>)> = if k == 0 {
            vec![
                (sais_ms, id),
                (sais_no_ms, id),
                (awry32, id),
                (sdsl_byte_32, id_skip),
            ]
        } else {
            vec![
                (sais_ms, min_no_remap),
                (sais_ms, min_remap),
                (sais_no_ms, min_no_remap),
                (sais_no_ms, min_remap),
                (awry32, min_no_remap),
                (awry32, min_remap),
                (sdsl_int_32, min_remap_skip),
                // TODO: SDSL without remap?
                // TODO: skip_zero_bytes for SDSL byte version
            ]
        };

        for (p, s) in params {
            tracing::warn!("Building UIndex with params {:?} {:?}", &*s, &*p);
            if let Some(u) = UIndex::<SV>::try_build_with_ranges(seq.clone(), &ranges, &*s, &*p) {
                let query_time = {
                    let _t = Timer::new("bench_positive").info();
                    u.bench_positive(&queries)
                };
                let mut stats = u.stats();
                stats.insert("sketch_params", Value::String(format!("{:?}", &*s)));
                stats.insert("index_params", Value::String(format!("{:?}", &*p)));
                stats.insert("query_length", Value::Number(Number::from(query_length)));
                stats.insert(
                    "query_time",
                    Value::Number(Number::from_f64(query_time).unwrap()),
                );
                all_stats.push(stats);
            }
        }

        // SIndex without sketching

        if k == 0 {
            tracing::warn!("Building plain-text SIndex");
            let u = SIndex::build(seq.clone(), 1, 1);

            let query_time = {
                let _t = Timer::new("bench_positive").info();
                u.bench_positive(&queries)
            };
            let mut stats = u.stats();
            stats.insert("index", Value::String("sparse SA".to_string()));
            stats.insert("index_params", Value::String("sparse SA".to_string()));
            stats.insert(
                "sketch_params",
                Value::String(format!("SIndex {{ k: {k}, l: {l} }}")),
            );
            stats.insert("sketch_l", Value::Number(Number::from(k)));
            stats.insert("sketch_k", Value::Number(Number::from(l)));
            stats.insert("query_length", Value::Number(Number::from(query_length)));
            stats.insert(
                "query_time",
                Value::Number(Number::from_f64(query_time).unwrap()),
            );
            all_stats.push(stats);
        }
        if k > 0 {
            tracing::warn!("Building minizer-space SIndex");
            let u = SIndex::build(seq.clone(), k, l);
            let query_time = {
                let _t = Timer::new("bench_positive").info();
                u.bench_positive(&queries)
            };
            let mut stats = u.stats();
            stats.insert("index", Value::String("SIndex".to_string()));
            stats.insert("index_params", Value::String("SIndex".to_string()));
            stats.insert(
                "sketch_params",
                Value::String(format!("SIndex {{ k: {k}, l: {l} }}")),
            );
            stats.insert("sketch_l", Value::Number(Number::from(k)));
            stats.insert("sketch_k", Value::Number(Number::from(l)));
            stats.insert("query_length", Value::Number(Number::from(query_length)));
            stats.insert(
                "query_time",
                Value::Number(Number::from_f64(query_time).unwrap()),
            );
            all_stats.push(stats);
        }
    }
}
