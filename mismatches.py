import os
import sys

l = 64
k = 8

total_mismatches = 0

with open("/home/blla008/Documents/u-index-rs/HG002") as file:
	for line in file:
		line_strip = line.rstrip()
			
		seqs = []
		for idx in range(0, len(line_strip), 128):
			seqs.append(line_strip[idx : idx + 128])
		
		f = open("/home/blla008/Documents/u-index-rs/sequence", "a")

		for seq in seqs:
			f.write(seq+'\n')
		f.close()
				
		os.system(r'cargo run -r --example bench -- --dna -k 8 -l 64  --text GRCh38.fna --patterns sequence')
		
		stats = open("/home/blla008/Documents/u-index-rs/stats.json", "r")	
		
		for line in stats:	
			if( '"query_mismatches"' in line ):
				x = line.find("query_mismatches")
				a = line.find(':', x)
				b = line.find(',', x)
				try:
					total_mismatches = total_mismatches + float( line[a+1:b] )
				except:
					total_mismatches = total_mismatches + float( line[a+1:b-1])
				
				print("Mismatches: ", str(total_mismatches) )
		
