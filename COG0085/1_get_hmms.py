import os
import sys
import subprocess
import re
import shlex
import pandas
import glob
import operator

input_list = open("/home/frankaylward/Desktop/marker_gene_benchmarking/genbank_curated_set/submodels/COG0085/rpob_model.list", "r")
pfam_folder = "/home/frankaylward/Desktop/marker_gene_benchmarking/hmm/Pfams/pfam_hmms/"
tigrfam_folder = "/home/frankaylward/Desktop/marker_gene_benchmarking/hmm/TIGRfams/TIGRFAMs_13.0_HMM/"
output_hmm1 = "/home/frankaylward/Desktop/marker_gene_benchmarking/genbank_curated_set/submodels/COG0085/pfam.HMM"
output_hmm2 = "/home/frankaylward/Desktop/marker_gene_benchmarking/genbank_curated_set/submodels/COG0085/tigrfam.HMM"

handle = open(output_hmm1, "w")
handle.close()
handle = open(output_hmm2, "w")
handle.close()

model_list = []
for i in input_list.readlines():
	model_list.append(i.rstrip())

for i in os.listdir(pfam_folder):
	model = re.sub(".hmm", "", i)
	if model in model_list:
		model_path = os.path.join(pfam_folder, i)
		cmd = "cat "+ model_path #+" > "+ output_hmm
		print cmd
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open(output_hmm1, "a"), stderr=open("err.txt", "w"))

for i in os.listdir(tigrfam_folder):
	model = re.sub(".HMM", "", i)
	if model in model_list:
		model_path = os.path.join(tigrfam_folder, i)
		cmd = "cat "+ model_path #+" > "+ output_hmm
		print cmd
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open(output_hmm2, "a"), stderr=open("err.txt", "w"))


output_hmm1_convert = re.sub(".HMM", ".hmm", output_hmm1)
cmd = "hmmconvert "+ output_hmm1
cmd2 = shlex.split(cmd)
subprocess.call(cmd2, stdout=open(output_hmm1_convert, "w"), stderr=open("err.txt", "w"))

output_hmm2_convert = re.sub(".HMM", ".hmm", output_hmm2)
cmd = "hmmconvert "+ output_hmm2
cmd2 = shlex.split(cmd)
subprocess.call(cmd2, stdout=open(output_hmm2_convert, "w"), stderr=open("err.txt", "w"))

cmd = "cat "+ output_hmm1_convert +" "+ output_hmm2_convert
cmd2 = shlex.split(cmd)
subprocess.call(cmd2, stdout=open("COG0085_submodels.hmm", "w"), stderr=open("err.txt", "w"))
