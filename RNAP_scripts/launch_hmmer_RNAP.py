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

input_folder = sys.argv[1]
output_dir = sys.argv[2]
output_file = sys.argv[3]
speci_db = "/home/frankaylward/Desktop/marker_gene_benchmarking/hmm/speci/RNAP.hmm"

# open and clear protein_statistics.txt file
#clear = open("protein_statistics.txt", "w")
#clear.close()
clear = open("error_file.txt", "w")
clear.close()

# define hmm launcher function
def hmm_launcher(folder):
	for files in os.listdir(input_folder):
		if files.endswith(".faa"):
			print files
			input_file = os.path.join(folder, files)	
			dom_output = re.sub(".faa", ".domout", files)
			speci_dom_output = os.path.join(output_dir, dom_output)

			# run against the RNAP models
			cmd = "hmmsearch --cpu 16 --domtblout "+ speci_dom_output +" "+ speci_db + " " + input_file
			print cmd
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdout=open("hmm.out", 'w'), stderr=open("error_file.txt", 'a'))

# end
hmm_launcher(input_folder)

################################################################
###### Loop through and parse the checkm HMM output ############
################################################################

def hmm_parser(folder, suffix, output):
	
	score_list = {}
	prot_list = []
	combined_output = open(output, "w")
	combined_output.write("protein\tacc\thit\tscore\tlength\tcategory\tspecies\tdomain\tphylum\torder\tp2\tp3\n")
	hits = []
	bit_dict = {}

	for filenames in os.listdir(folder):
		if filenames.endswith(suffix):

			acc = re.sub(suffix, "", filenames)

			f = open(folder+"/"+filenames, 'r')
			o = open(folder+"/"+filenames+".parsed", 'w')
			hit_dict = {}
			bit_dict = defaultdict(int)
			hit_type = {}
			marker_dict = {}

			for line in f.readlines():
				if line.startswith("#"):
					pass
				else:
					newline = re.sub( '\s+', '\t', line)
					list1 = newline.split('\t')
					ids = list1[0]
					hit = list1[2]
					acc2 = list1[3]
					if "COG" in hit:
						pass
					else:
						hit = acc2

					bit_score = list1[7]
					score = float(bit_score)

					if score > 500:
						if score > bit_dict[ids]:
							hit_dict[ids] = hit
							bit_dict[ids] = score
					else:
						print ids, hit, score

			bit_sorted = sorted(bit_dict.items(), key=operator.itemgetter(1), reverse=True)
			output_list = []
			for item in bit_sorted:
				entry = item[0]
				output_list.append(entry +"\t"+ str(hit_dict[entry]) +"\t"+ str(bit_dict[entry]) )

			hit_profile = defaultdict(int)
			done = []
			for line in output_list:
				line1 = line.rstrip()
				tabs = line1.split("\t")
				ids = tabs[0]
				hits.append(ids)
				cog = tabs[1]
				score = tabs[2]
				nr = acc +"_"+ cog

				if nr in done:
					pass
					#combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\tNH\n")
					#o.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\tNH\n")
				else:
					combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\tBH\n")
					o.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\tBH\n")
					done.append(nr)
			o.close()

#parse speci outputs
speci_df = hmm_parser(output_dir, ".domout", output_file)




