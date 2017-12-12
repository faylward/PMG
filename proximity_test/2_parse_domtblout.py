import os
import sys
import subprocess
import re
import shlex
import pandas
import glob
import operator
from collections import defaultdict
from Bio import SeqIO

output_dir = "/home/frankaylward/Desktop/marker_gene_benchmarking/metagenomes/metagenomes_5KB/parks_et_al_raw/"
#genomes = open("/home/frankaylward/Desktop/marker_gene_benchmarking/genbank_curated_set/good_quality_genomes.list", "r")

################################################################
###### Loop through and parse the checkm HMM output ############
################################################################

def hmm_parser(folder, suffix, output):
	
	score_list = {}
	prot_list = []
	hits = []
	bit_dict = {}
	len_dict = {}
	
	for filenames in os.listdir(folder):
		if filenames.endswith(suffix):

			acc = re.sub(suffix, "", filenames)
			print acc
			f = open(folder+"/"+filenames, 'r')
			o = open(folder+"/"+filenames+".parsed", 'w')
			hit_dict = {}
			bit_dict = defaultdict(int)

			start_dict = {}
			end_dict = {}

			hit_type = {}
			marker_dict = {}

			for line in f.readlines():
				if line.startswith("#"):
					pass
				else:
					newline = re.sub( '\s+', '\t', line)
					list1 = newline.split('\t')
					ids = list1[0]
					tlen = list1[2]
					hit = list1[3]
					bit_score = list1[7]
					score = float(bit_score)

					id_hit = ids +"|"+ hit
					len_dict[id_hit] = tlen
					start = int(list1[15])
					end = int(list1[16])
						
					if id_hit in start_dict:
						if start_dict[id_hit] > start:
							start_dict[id_hit] = start
						if end_dict[id_hit] < end:
							end_dict[id_hit] = end
					else:
						start_dict[id_hit] = start
						end_dict[id_hit] = end

					if score > 0:

						if score > bit_dict[id_hit]:
							hit_dict[id_hit] = hit
							bit_dict[id_hit] = score
					else:
						print ids, hit, score

			bit_sorted = sorted(bit_dict.items(), key=operator.itemgetter(1), reverse=True)
			for item in bit_sorted:
				entry = item[0]
				entries = entry.split("|")
				item = entries[0]
				o.write(item +"\t"+ str(hit_dict[entry]) +"\t"+ len_dict[entry] +"\t"+ str(bit_dict[entry]) +"\t"+ str(start_dict[entry]) +"\t"+ str(end_dict[entry]) +"\n")

			o.close()
	return entry

#parse speci outputs
speci_df = hmm_parser(output_dir, ".hmm_dom_out", "all_hmm_domout.speci.txt")




