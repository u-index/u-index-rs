#![feature(lazy_get)]
#![allow(unused)]
use std::{any::type_name, cell::Cell, collections::HashMap, sync::LazyLock};

use clap::Parser;
use packed_seq::{AsciiSeqVec, PackedSeqVec, SeqVec};
use sdsl_lite_fm::*;
use serde_json::{Number, Value};
use tracing::info;
use uindex::{
    bench::gen_query_positions,
    indices::{DivSufSortSa, FmAwryParams, FmSdslParams, LibSaisSa},
    s_index::SIndex,
    sketchers::{IdentityParams, MinimizerParams},
    utils::{read_chromosomes, Timer, INIT_TRACE},
    IndexBuilder, SketcherBuilder, UIndex,
};

#[derive(Parser)]
struct Args {
    /// Run experiments in subprocesses.
    #[clap(long)]
    ext: bool,
    #[clap(long)]
    run_id: Option<usize>,
}

static ARGS: LazyLock<Args> = LazyLock::new(|| Args::parse());

thread_local! {
    static ID: Cell<usize> = Cell::new(0);
}

fn main() {
    color_backtrace::install();
    *INIT_TRACE;

    if ARGS.run_id.is_none_or(|run_id| run_id < 1000) {
        let (seq, _ranges) = read_chromosomes::<PackedSeqVec>(1);

        let query_length = 512;
        let num_queries = 10000;
        let queries = gen_query_positions(seq.len(), query_length, num_queries);

        let mut all_stats = vec![];
        run::<PackedSeqVec>(
            &mut all_stats,
            seq,
            query_length,
            &queries,
            &[(0, 0), (4, 32), (8, 64), (16, 128), (28, 256)],
        );

        // Write all_stats.
        if ARGS.run_id.is_none() {
            let stats_string = serde_json::to_string(&all_stats).unwrap();
            std::fs::write("stats.json", stats_string).unwrap();
        }
    }

    ID.set(1000);
    if ARGS
        .run_id
        .is_none_or(|run_id| 1000 <= run_id && run_id < 2000)
    {
        let seq = std::fs::read("english.200MB").unwrap();

        let query_length = 128;
        let num_queries = 10000;
        let queries = gen_query_positions(seq.len(), query_length, num_queries);

        let mut all_stats = vec![];
        run::<Vec<u8>>(
            &mut all_stats,
            seq,
            query_length,
            &queries,
            &[(0, 0), (2, 20), (3, 24), (4, 32), (8, 64)],
        );

        // Write all_stats.
        if ARGS.run_id.is_none() {
            let stats_string = serde_json::to_string(&all_stats).unwrap();
            std::fs::write("stats-english.json", stats_string).unwrap();
        }
    }

    ID.set(2000);
    if ARGS
        .run_id
        .is_none_or(|run_id| 2000 <= run_id && run_id < 3000)
    {
        let seq = std::fs::read("proteins.200MB").unwrap();

        let query_length = 128;
        let num_queries = 10000;
        let queries = gen_query_positions(seq.len(), query_length, num_queries);

        let mut all_stats = vec![];
        run::<Vec<u8>>(
            &mut all_stats,
            seq,
            query_length,
            &queries,
            &[(0, 0), (2, 20), (3, 24), (4, 32), (8, 64)],
        );

        // Write all_stats.
        if ARGS.run_id.is_none() {
            let stats_string = serde_json::to_string(&all_stats).unwrap();
            std::fs::write("stats-proteins.json", stats_string).unwrap();
        }
    }

    if let Some(run_id) = ARGS.run_id {
        panic!("Did non find run id {run_id}");
    }
}

fn run<SV: SeqVec>(
    all_stats: &mut Vec<HashMap<&str, Value>>,
    seq: SV,
    query_length: usize,
    queries: &Vec<(usize, usize)>,
    kl: &[(usize, usize)],
) {
    tracing::warn!("{}", type_name::<SV>());
    let ranges = vec![0..seq.len()];

    for &(k, l) in kl {
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
                (sdsl_byte_32, id),
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
            let stats = run_fn(|| {
                tracing::warn!("Building UIndex with params {:?} {:?}", &*s, &*p);

                let u = UIndex::<SV>::try_build_with_ranges(seq.clone(), &ranges, &*s, &*p)?;
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
                Some(stats)
            });
            if let Some(stats) = stats {
                all_stats.push(stats);
            }
        }

        // SIndex without sketching

        if k == 0 {
            let stats = run_fn(|| {
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
                stats.insert("sketch_l", Value::Number(Number::from(l)));
                stats.insert("sketch_k", Value::Number(Number::from(k)));
                stats.insert("query_length", Value::Number(Number::from(query_length)));
                stats.insert(
                    "query_time",
                    Value::Number(Number::from_f64(query_time).unwrap()),
                );
                Some(stats)
            });
            if let Some(stats) = stats {
                all_stats.push(stats);
            }
        }
        if k > 0 {
            let stats = run_fn(|| {
                tracing::warn!("Building minizer-space SIndex");
                let u = SIndex::build(seq.clone(), k, l);
                let query_time = {
                    let _t = Timer::new("bench_positive").info();
                    u.bench_positive(&queries)
                };
                let mut stats = u.stats();
                stats.insert("index", Value::String("sparse SA".to_string()));
                stats.insert("index_params", Value::String("SIndex".to_string()));
                stats.insert(
                    "sketch_params",
                    Value::String(format!("SIndex {{ k: {k}, l: {l} }}")),
                );
                stats.insert("sketch_l", Value::Number(Number::from(l)));
                stats.insert("sketch_k", Value::Number(Number::from(k)));
                stats.insert("query_length", Value::Number(Number::from(query_length)));
                stats.insert(
                    "query_time",
                    Value::Number(Number::from_f64(query_time).unwrap()),
                );
                Some(stats)
            });
            if let Some(stats) = stats {
                all_stats.push(stats);
            }
        }
    }
}

fn run_fn(
    f: impl Fn() -> Option<HashMap<&'static str, Value>>,
) -> Option<HashMap<&'static str, Value>> {
    if ARGS.ext {
        ID.set(ID.get() + 1);
        if let Some(run_id) = ARGS.run_id {
            if run_id == ID.get() {
                info!("RUNNING FOR ID {}", ID.get());
                let mut result = f();
                if let Some(result) = result.as_mut() {
                    result.insert("id", Value::Number(Number::from(ID.get())));
                }
                println!("{}", serde_json::to_string(&result).unwrap());
                info!("RUNNING FOR ID {} DONE", ID.get());
                info!("OUTPUT {}", serde_json::to_string(&result).unwrap());
                std::process::exit(0);
            } else {
                // This is not the run.
                None
            }
        } else {
            // Launch subprocess to run this experiment.
            let program = std::env::current_exe().unwrap();
            info!("LAUNCHING subprocess for ID {} {:?}", ID.get(), program);
            let output = std::process::Command::new(&program)
                .arg("--ext")
                .arg("--run-id")
                .arg(&format!("{}", ID.get()))
                .stderr(std::process::Stdio::inherit())
                .output()
                .unwrap();
            info!("READING   subprocess for ID {} {:?}", ID.get(), program);
            // FIXME: Just leaking all the output so that we can get static strings to it.
            let output = Box::new(output.stdout).leak();
            // info!("READING   OUTPUT {}", std::str::from_utf8(output).unwrap());
            let result = serde_json::from_slice(output).unwrap();
            result
        }
    } else {
        // Just run directly here.
        f()
    }
}
