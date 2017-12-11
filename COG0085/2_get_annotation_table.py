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

hmm = open("/home/frankaylward/Desktop/marker_gene_benchmarking/genbank_curated_set/submodels/COG0085/COG0085_vs_submodels.hmmout", "r")
marker_out = open("/home/frankaylward/Desktop/marker_gene_benchmarking/genbank_curated_set/merged_COG0085_annote.txt", "r")
full_dict = defaultdict(lambda: defaultdict(float))
done = {}
done1 = {}
protein_list = []
features = ['qlen']

df = pandas.DataFrame()
for i in hmm.readlines():
	if i.startswith("#"):
		pass
	else:
		line = i.rstrip()
		line2 = re.sub("\s+", "\t", line)	
		#line3 = re.sub("\s", "\t", line)
		tabs = line2.split("\t")
		#print tabs
		protein = tabs[0]
		hit = tabs[4]
		qlen = tabs[2]
		score = tabs[7]
		full = protein +"_"+ hit

		if full in done1:
			pass
		else:
			protein_list.append(protein)
			done1[full] = full
			full_dict[hit][protein] = score
			#print qlen
			features.append(hit)
			
for i in marker_out.readlines():
	line = i.rstrip()
	tabs = line.split("\t")
	hit = tabs[5]
	protein = tabs[0]
	if hit == "COG0085": #and protein in protein_list:
		hit_type = tabs[13]
		score = tabs[8]
		domain = tabs[3]
		phylum = tabs[4]
		length = tabs[7]

		full_dict["category"][protein] = hit_type
		full_dict["COG0085_score"][protein] = score		
		full_dict["domain"][protein] = domain
		full_dict["phylum"][protein] = phylum
		full_dict["qlen"][protein] = length
	else:
		print protein


for j in full_dict:
	s1 = pandas.DataFrame(pandas.Series(full_dict[j], name = j))
	df = pandas.concat([df, s1], axis=1)

print len(set(protein_list))

feature_set = list(set(features))
column_names = ['domain', 'phylum', 'category', 'COG0085_score'] + feature_set
df = df[column_names]
name2 = 'protein_annotations.txt'
df.fillna(0, inplace=True, axis=1)
df.to_csv(name2, sep='\t')
		
	
