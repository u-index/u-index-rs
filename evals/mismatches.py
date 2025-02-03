#!/usr/bin/env python3
from pathlib import Path
import os
import sys

block = 128  # change this to be the same size as l

mismatches = 0
all_stats = ""
with open("hg-long-reads") as file:
    for _ in range(10):
        line = file.readline()
        l = line.rstrip()
        seqs = []
        for idx in range(0, len(l), block):
            seqs.append(l[idx : idx + block])

        f = open("/tmp/patterns", "w")

        for seq in seqs:
            f.write(seq + "\n")
        f.close()

        os.system(
            r"cargo run -r --example bench -- --dna -k 16 -l 128 --text human-genome.fa --patterns /tmp/patterns"
        )

        stats = Path("stats.json").read_text()
        all_stats += stats + "\n"
        print(stats)

    #     for line in stats:
    #         if '"query_matches":0.0' in line:
    #             mismatches = mismatches + 1

    # mis = open("output", "a")
    # mis.write(mismatches)

print(all_stats)

Path("all_stats.json").write_text(all_stats)
