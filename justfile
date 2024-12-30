py:
    source .env/bin/activate && maturin develop -r

py-init:
    python3 -m venv .env
    ln -sfn ../.env/lib/python3.13/site-packages/uindex py/uindex
    ln -sfn ../human-genome.fa py/human-genome.fa

perf:
    cargo build -r --example bench
    perf record target/release/examples/bench
    perf report

flame:
    cargo flamegraph --open --example bench
