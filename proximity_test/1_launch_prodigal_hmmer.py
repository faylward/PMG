import argparse
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

input_folder = "/home/frankaylward/Desktop/marker_gene_benchmarking/metagenomes/metagenomes_5KB/parks_et_al_raw/"
db = "/home/frankaylward/Desktop/marker_gene_benchmarking/hmm/speci/all.hmm"

#clear = open("protein_statistics.txt", "w")
#clear.close()

#final_set = ['ERX552293_contigs.10KB.fna', 'SRX290715_contigs.10KB.fna', 'SRX894049_contigs.10KB.fna']
#final_set = ['SRX035182_contigs.10KB.fna', 'SRX516605_contigs.10KB.fna']

# define hmm launcher function
for folders in os.listdir(input_folder):
	if folders.endswith(".fna"):

		#if folders in final_set:

		print folders
		input_file = os.path.join(input_folder, folders)
		protein_file = re.sub(".fna", ".faa", input_file)
		gff_file = re.sub(".fna", ".gff", input_file)

		cmd = "prodigal -i "+ input_file +" -p meta -f gff -o "+ gff_file +" -a "+ protein_file 
		#print cmd
		cmd2 = shlex.split(cmd)
		#subprocess.call(cmd2, stdout=open("prodigal.out", 'w'), stderr=open("prodigal.err", 'a'))

		hmm_out = re.sub(".fna", ".hmm_out", input_file)
		dom_out = re.sub(".fna", ".hmm_dom_out", input_file)

		cmd = cmd = "hmmsearch --cpu 16 --tblout " + hmm_out + " --domtblout "+ dom_out + " " + db + " " + protein_file
		#print cmd
		cmd2 = shlex.split(cmd)
		#subprocess.call(cmd2, stdout=open("hmmer.out", 'w'), stderr=open("hmmer.err", 'a'))

		#for i in SeqIO.parse(protein_file, "r")

		# set up seqtk command and call it
		cmd = "seqtk comp " + protein_file
		cmd2 = shlex.split(cmd)
		#subprocess.call(cmd2, stdout=open("protein_statistics.txt", 'a'), stderr=open("error_file.txt", 'a'))


protein_lengths = {}
# get dictionary of protein lengths
protein_stats = open("protein_statistics.txt", "r")
for i in protein_stats.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	name = tabs[0]
	length = tabs[1]
	protein_lengths[name] = length

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
	df = pandas.DataFrame()

	for filenames in os.listdir(folder):
		if filenames.endswith(suffix):

			acc = re.sub(suffix, "", filenames)
			if 1 == 1:

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

						bit_score = list1[5]
						score = float(bit_score)

						if score > 0:

							if score > bit_dict[ids]:
								hit_dict[ids] = hit
								bit_dict[ids] = score
						else:
							print ids, hit, score

				bit_sorted = sorted(bit_dict.items(), key=operator.itemgetter(1), reverse=True)
				output_list = []
				for item in bit_sorted:
					entry = item[0]
					output_list.append(entry +"\t"+ str(hit_dict[entry]) +"\t"+ str(bit_dict[entry]) +"\t"+ str(protein_lengths[entry]))
					#o.write(entry +"\t"+ str(hit_dict[entry]) +"\t"+ str(bit_dict[entry]) +"\t"+ str(protein_lengths[entry]) +"\n")
				#o.close()

				#parsed = open(folder+"/"+filenames+".parsed", 'r')
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
					#acc2 = taxon_dict2[acc]

					if nr in done:
						combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\t"+ str(protein_lengths[ids]) +"\tNH\t"+ "\n")#taxon_dict[acc] +"\n")
						o.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\t"+ str(protein_lengths[ids]) +"\tNH\t"+ "\n") #taxon_dict[acc] +"\n")
					else:
						combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\t"+ str(protein_lengths[ids]) +"\tBH\t"+ "\n")#taxon_dict[acc] +"\n")
						o.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ score +"\t"+ str(protein_lengths[ids]) +"\tBH\t"+ "\n")#taxon_dict[acc] +"\n")
						done.append(nr)

				o.close()

				s1 = pandas.DataFrame(pandas.Series(hit_profile, name = acc))
				df = pandas.concat([df, s1], axis=1)
	return df

#parse speci outputs
speci_df = hmm_parser(input_folder, ".hmm_out", "all_hmm_out.speci.txt")
name2 = 'hmm_profile.speci.txt'
speci_df.fillna(0, inplace=True, axis=1)
df2 = speci_df.transpose()
df2.to_csv(name2, sep='\t')


