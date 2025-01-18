use packed_seq::{PackedSeqVec, SeqVec};
use sdsl_lite_fm::*;
use serde_json::{Number, Value};
use uindex::{
    bench::gen_query_positions,
    indices::{DivSufSortSa, FmAwryParams, FmSdslParams, LibSaisSa},
    s_index::SIndex,
    sketchers::{IdentityParams, MinimizerParams},
    utils::{read_chromosomes, Timer},
    IndexBuilder, SketcherBuilder, UIndex,
};

fn main() {
    color_backtrace::install();
    *INIT_TRACE;

    let mut all_stats = vec![];
    let query_length = 512;
    let num_queries = 10000;
    let queries = gen_query_positions(seq.len(), query_length, num_queries);

    let (seq, ranges) = read_chromosomes::<PackedSeqVec>(1);

    for (k, l) in [
        (4, 32),
        // (8, 64), (15, 128), (28, 256)
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

        let params: Vec<(&dyn IndexBuilder, &dyn SketcherBuilder)> = vec![
            // (sais_ms, id),
            (sais_ms, min_no_remap),
            // (sais_ms, min_remap),
            // (sais_no_ms, id),
            // (sais_no_ms, min_no_remap),
            // (sais_no_ms, min_remap),
            // (awry32, id),
            // (awry32, min_no_remap),
            // (awry32, min_remap),
            // (sdsl_byte_32, id_skip),
            // (sdsl_int_32, min_remap_skip),
            // // TODO: SDSL without remap?
            // // TODO: skip_zero_bytes for SDSL byte version
        ];

        for (p, s) in params {
            if let Some(u) = UIndex::try_build_with_ranges(seq.clone(), &ranges, &*s, &*p) {
            tracing::warn!("Building UIndex with params {:?} {:?}", &*s, &*p);
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

        if false {
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
        {
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

    // Write all_stats.
    let stats_string = serde_json::to_string(&all_stats).unwrap();
    std::fs::write("stats.json", stats_string).unwrap();
}
