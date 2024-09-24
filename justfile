py:
    source .env/bin/activate && maturin develop -r

py-init:
    python3 -m venv .env
    ln -sfn ../.env/lib/python3.12/site-packages/uindex py/uindex
    ln -sfn ../human-genome.fa py/human-genome.fa
