import os, sys, re, operator, subprocess, shlex, glob
from collections import defaultdict

input_dir  = sys.argv[1]
output_dir = sys.argv[2]

folder_list = []
for folders in os.listdir(input_dir):
	subfolders = os.path.join(input_dir, folders)
	#print subfolders
	file_list = os.listdir(subfolders)
	list1 = [x for x in file_list if "protein.faa.gz" in x]
	list2 = [x for x in file_list if ".faa" in x]
	if len(list1) < 1:
		# If this is happening, then no protein file could be detected, and we have to predict our own.
		list3 = [x for x in file_list if "genomic.fna.gz" in x or "genomic.fna" in x]
		genome_file = list3[0]
		genome_path = os.path.join(subfolders, genome_file)

		if genome_path.endswith(".gz"):
			# unzip genome file
			cmd = "gunzip "+ genome_path
			cmd1 = shlex.split(cmd)
			subprocess.call(cmd1, stdout=open("stdout.txt", "w"), stderr=open("stderr.txt", "w"))
			genome_path = re.sub(".gz", "", genome_path)

		# create names and paths for output files
		protein_output = os.path.join(output_dir, re.sub(".fna", ".prodigal.faa", genome_file))
		gff_output = os.path.join(output_dir, re.sub(".fna", ".prodigal.gff", genome_file))

		# Run prodigal to predict genes/proteins
		cmd = "prodigal -f gff -i "+ genome_path +" -a "+ protein_output +" -o "+ gff_output
		cmd1 = shlex.split(cmd)
		subprocess.call(cmd1, stdout=open("prodigal.out", "w"), stderr=open("prodigal.err", "w"))

	else:
		# if this happens, then a protein file could be identified, and we need to move it to the final output folder
		list3 = [x for x in file_list if "_protein.faa.gz" in x or "protein.faa" in x]
		list4 = [x for x in list3 if not "_from" in x]

		genome_file = list4[0]
		genome_path = os.path.join(subfolders, genome_file)

		#if genome_path.endswith(".gz"):
		#	# unzip genome file
		#	cmd = "gunzip "+ genome_path
		#	cmd1 = shlex.split(cmd)
		#	subprocess.call(cmd1, stdout=open("stdout.txt", "w"), stderr=open("stderr.txt", "w"))
		#	genome_path = re.sub(".gz", "", genome_path)

		gff_list = [x for x in file_list if ".gff" in x]
		gff_file = gff_list[0]
		gff_path = os.path.join(subfolders, gff_file)

		cmd = "cp "+ genome_path +" "+ output_dir
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open("prodigal.out", "w"), stderr=open("prodigal.err", "w"))	

		cmd = "cp "+ gff_path +" "+ output_dir
		cmd2 = shlex.split(cmd)
		subprocess.call(cmd2, stdout=open("prodigal.out", "w"), stderr=open("prodigal.err", "w"))



