import os
import sys

block = 128 # change this to be the same size as l

mismatches = 0
with open("/home/blla008/Documents/u-index-rs/HG002") as file:
	while file:
		line = file.readline()
		l = line.rstrip()
		
		seqs = []
		for idx in range(0, len(l), block):
			seqs.append(l[idx : idx + block])
				
		f = open("/home/blla008/Documents/u-index-rs/sequence", "w")

		for seq in seqs:
			f.write(seq+'\n')
		f.close()
			
		os.system(r'cargo run -r --example bench -- --dna -k 16 -l 128 --text GRCh38.fna --patterns sequence')
	
		stats = open("/home/blla008/Documents/u-index-rs/stats.json", "r")
		
		for line in stats:
			if ('"query_matches":0.0' in line ):
				mismatches = mismatches + 1
		
		
	mis = open("/home/blla008/Documents/u-index-rs/16,128", 'a')
	mis.write( mismatches ")
