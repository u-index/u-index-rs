#!/usr/bin/env python3
import os
import sys

block = 128  # change this to be the same size as l

mismatches = 0
with open("long-reads") as file:
    for _ in range(100):
        line = file.readline()
        l = line.rstrip()
        seqs = []
        for idx in range(0, len(l), block):
            seqs.append(l[idx : idx + block])

        f = open("patterns", "w")

        for seq in seqs:
            f.write(seq + "\n")
        f.close()

        os.system(
            r"cargo run -r --example bench -- --dna -k 16 -l 128 --text human-genome.fa --patterns patterns"
        )

        stats = open("stats.json", "r")
        print(stats)

        for line in stats:
            if '"query_matches":0.0' in line:
                mismatches = mismatches + 1
        break

    mis = open("output", "a")
    mis.write(mismatches)
