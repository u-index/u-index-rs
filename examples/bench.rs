#![feature(lazy_get)]
#![allow(unused)]
use std::{any::type_name, cell::Cell, collections::HashMap, path::PathBuf, sync::LazyLock};

use clap::Parser;
use packed_seq::{AsciiSeqVec, PackedSeqVec, SeqVec};
use sdsl_lite_fm::*;
use serde_json::{Number, Value};
use tracing::info;
use uindex::{
    bench::gen_queries,
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

    /// Plain-text input file.
    #[clap(long)]
    text: Option<PathBuf>,
    /// Query patterns. One per line.
    #[clap(long)]
    patterns: Option<PathBuf>,
    /// Output json file with statistics.
    #[clap(long)]
    output: Option<PathBuf>,

    #[clap(short)]
    k: Option<usize>,
    #[clap(short)]
    l: Option<usize>,
}

static ARGS: LazyLock<Args> = LazyLock::new(|| Args::parse());

thread_local! {
    static ID: Cell<usize> = Cell::new(0);
}

fn main() {
    color_backtrace::install();
    *INIT_TRACE;

    let mut kls = &[(0, 0), (4, 32), (8, 64), (16, 128), (28, 256)][..];

    if ARGS.k.is_some() ^ ARGS.l.is_some() {
        panic!("Both or none of k and l must be specified");
    }

    let kl = [(ARGS.k.unwrap_or_default(), ARGS.l.unwrap_or_default())];
    if ARGS.k.is_some() {
        kls = &kl[..];
    }

    if ARGS.patterns.is_some() ^ ARGS.text.is_some() {
        panic!("Both or none of --text and --patterns must be specified.");
    }

    if let Some(text) = &ARGS.text {
        let seq = std::fs::read(text).unwrap();

        let pattern_data = std::fs::read(ARGS.patterns.as_ref().unwrap()).unwrap();
        let queries = pattern_data
            .split(|&c| c == b'\n')
            .map(|x| x.to_vec())
            .collect::<Vec<_>>();
        let query_length = 0;

        let mut all_stats = vec![];
        run::<Vec<u8>>(&mut all_stats, &seq, query_length, &queries, kls);

        let stats_string = serde_json::to_string(&all_stats).unwrap();
        let path = PathBuf::from("stats.json");
        let output = ARGS.output.as_ref().unwrap_or(&path);
        std::fs::write(output, stats_string).unwrap();

        return;
    }

    if ARGS.run_id.is_none_or(|run_id| run_id < 1000) {
        let (seq, _ranges) = read_chromosomes::<PackedSeqVec>(1);

        let query_length = 512;
        let num_queries = 10000;
        let queries = gen_queries::<PackedSeqVec>(&seq, query_length, num_queries);

        let mut all_stats = vec![];
        run::<PackedSeqVec>(&mut all_stats, &seq, query_length, &queries, kls);

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
        let queries = gen_queries::<Vec<u8>>(&seq, query_length, num_queries);

        let mut all_stats = vec![];
        run::<Vec<u8>>(&mut all_stats, &seq, query_length, &queries, kls);

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
        let queries = gen_queries::<Vec<u8>>(&seq, query_length, num_queries);

        let mut all_stats = vec![];
        run::<Vec<u8>>(&mut all_stats, &seq, query_length, &queries, kls);

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

fn run<'s, SV: SeqVec>(
    all_stats: &mut Vec<HashMap<&str, Value>>,
    seq: &'s SV,
    query_length: usize,
    queries: &[SV],
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

                let rss0 = max_rss();

                let u = UIndex::<SV>::try_build_with_ranges(&seq, &ranges, &*s, &*p)?;
                let rss1 = max_rss();
                let query_time = {
                    let _t = Timer::new("bench_positive").info();
                    u.bench(queries)
                };
                let rss2 = max_rss();
                let mut stats = u.stats();
                stats.insert("rss0", Value::Number(Number::from(rss0)));
                stats.insert("rss1", Value::Number(Number::from(rss1)));
                stats.insert("rss2", Value::Number(Number::from(rss2)));
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
                let rss0 = max_rss();
                let u = SIndex::build(seq, 1, 1);
                let rss1 = max_rss();

                let query_time = {
                    let _t = Timer::new("bench_positive").info();
                    u.bench(queries)
                };
                let rss2 = max_rss();
                let mut stats = u.stats();
                stats.insert("rss0", Value::Number(Number::from(rss0)));
                stats.insert("rss1", Value::Number(Number::from(rss1)));
                stats.insert("rss2", Value::Number(Number::from(rss2)));
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
                let rss0 = max_rss();
                let u = SIndex::build(seq, k, l);
                let rss1 = max_rss();
                let query_time = {
                    let _t = Timer::new("bench_positive").info();
                    u.bench(&queries)
                };
                let rss2 = max_rss();
                let mut stats = u.stats();
                stats.insert("rss0", Value::Number(Number::from(rss0)));
                stats.insert("rss1", Value::Number(Number::from(rss1)));
                stats.insert("rss2", Value::Number(Number::from(rss2)));
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

fn max_rss() -> usize {
    let rusage = unsafe {
        let mut rusage = std::mem::MaybeUninit::uninit();
        libc::getrusage(libc::RUSAGE_SELF, rusage.as_mut_ptr());
        rusage.assume_init()
    };
    // On linux, the returned value is in kB.
    rusage.ru_maxrss as usize * 1024
}
