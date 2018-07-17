import os
import sys
import subprocess
import re
import shlex
import pandas
import glob
import operator
import numpy as np
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

working_dir = sys.argv[1]
info = open("complete_genome_reps.txt", "r")

merged = open("RNAP.txt", "w")
merged.write("protein\tacc\tspecies\tdomain\tphylum\tfamily\tgenus\thit\tcategory\tlength\tscore\talign_length\tnum_hits\tall_proteins\talignment_locations\n")

merged_proteins = open("RNAP.faa", "w")
final_proteins = []

# how many genes to look on either side of the initial RNAP hits?
prox = int(50)

# get taxonomic info
taxon_dict = {}
for j in info:
	line = j.rstrip("\n")
	tabs = line.split("\t")
	ftp = tabs[19]
	slash = ftp.split("/")
	name = slash[9]
	tax = [x for i,x in enumerate(tabs) if i in [7, 22, 23, 24, 25]]

	taxonomic = "\t".join(tax)
	taxon_dict[name] = taxonomic


marker_tally = defaultdict(int)
exceptions = open("exceptions.txt", "w")
tally = 0

#######################################################
########## Define function to parse GFF files #########
#######################################################
def gff_parser(path_to_gff_file):
	region2proteins = defaultdict(list)
	protein2region = {}
	gff_handle = open(path_to_gff_file, "r")

	if "prodigal" in path_to_gff_file:
		for n in gff_handle.readlines():
			if n.startswith("#"):
				pass
			else:
				line = n.rstrip()
				tabs = line.split("\t")
				if tabs[2] == "CDS":
					region = tabs[0]
					start = tabs[3]
					end = tabs[4]
					full = tabs[8]
					semi = full.split(";")
					name = semi[0]
					underscore = name.split("_")
					num = underscore[1]
					prot = region +"_"+ str(num)

					region2proteins[region].append(prot)
					protein2region[prot] = region
	else:
		for n in gff_handle.readlines():
			if n.startswith("#"):
				pass
			else:
				line = n.rstrip()
				tabs = line.split("\t")
				if tabs[2] == "CDS":
					region = tabs[0]
					start = tabs[3]
					end = tabs[4]
					full = tabs[8]
					if "pseudo=true" in full:
						pass
					else:
						semi = full.split(";")
						for m in semi:
							if "protein_id" in m:
								prot = re.sub("protein_id=", "", m)

						region2proteins[region].append(prot)
						protein2region[prot] = region

	gff_handle.close()
	return region2proteins, protein2region


################################################################
########## Define function for parsing HMMER3 output ###########
################################################################
def parse_domout(path_to_parsed_hmmfile, acc, protein_dict):
	parsed = open(path_to_parsed_hmmfile, "r")
	done = {}
	protein2coords = defaultdict(list)
	protein2align_length = {}

	main_hit = "NAN"
	rnap_hits = []
	protein2cog = defaultdict(lambda:"NA")
	protein2acc = {}
	protein2score = {}
	protein2category = {}
	protein2length = {}

	for n in parsed.readlines():
		line = n.rstrip()
		tabs = line.split("\t")
		protein = tabs[0]
		annot = tabs[2]
		if annot == "COG0085":
			rnap_hits.append(protein)
			id_hit = protein +"|"+ annot
			hmm_score = tabs[5]
			category = tabs[6]

			start = int(tabs[3])
			end =   int(tabs[4])

			record = protein_dict[protein]
			prot_length = len(record.seq)

			nr = acc +"_"+ annot

			protein2cog[protein]    = annot
			protein2acc[protein]    = acc
			protein2score[protein]  = hmm_score
			protein2length[protein] = prot_length

			protein2coords[id_hit].append(start)
			protein2coords[id_hit].append(end)

			align_length = abs(end - start)
			protein2align_length[id_hit] = align_length

			if category in "BH":
				main_hit = protein
				protein2dups[id_hit]

			protein2category[protein] = category

	parsed.close()
	return main_hit, protein2cog, protein2acc, protein2score, protein2length, protein2category, protein2coords, protein2align_length


for i in os.listdir(working_dir):
	if i.endswith(".faa"):
		#print i
		protein2dups = defaultdict(lambda:"single_besthit")
		if "prodigal" in i:
			protein_file = os.path.join(working_dir, i)
			gff_file = re.sub(".faa", ".gff", protein_file)
			domout = re.sub(".faa", ".domout", protein_file)
			parsed = re.sub(".faa", ".domout.parsed", protein_file)
			acc = re.sub("_genomic.prodigal.faa", "", i)	
		else:
			protein_file = os.path.join(working_dir, i)
			gff_file = re.sub("_protein.faa", "_genomic.gff", protein_file)
			domout = re.sub("_protein.faa", "_protein.domout", protein_file)
			parsed = re.sub("_protein.faa", "_protein.domout.parsed", protein_file)
			acc = re.sub("_protein.faa", "", i)	

		# parse gff file and get coordinate dictionaries
		region2proteins, protein2region = gff_parser(gff_file)

		# get a dictionary of protein sequences
		seq_handle = open(protein_file, "r")
		seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))

		# parse domout file and get protein hits and coordinates
		rnap, protein2cog, protein2acc, protein2score, protein2length, protein2category, protein2coords, protein2align_length = parse_domout(parsed, acc, seq_dict)

		if rnap == "NAN":
			print acc
		else:

			already_done = []
			num_proteins = defaultdict(lambda:int(1))
			prot2protlist = defaultdict(list)
			prot2loc = defaultdict(list)

			prot2protlist[rnap].append(rnap)

			id_hit1 = rnap +"|COG0085"
			range1 = protein2coords[id_hit1]

			r1 = range(range1[0], range1[1])
			meanloc1 = np.mean(range1)
			prot2loc[rnap].append(meanloc1)

			region = protein2region[rnap]

			protein_list = region2proteins[region]
			index = int(protein_list.index(rnap))
			length = len(protein_list)

			if index + prox > length:
				right_bound = length
			else:
				right_bound = index + prox

			if index - prox < 0:
				left_bound = 0
			else:
				left_bound = index - prox

			orf_set = set(protein_list[left_bound:right_bound])

			if rnap in orf_set:
				orf_set.remove(rnap)
			else:
				exceptions.write(start +"\t"+ end +"\t"+ left_bound +"\t"+ right_bound +"\t"+ region_length +"\n")

			for d in orf_set:
				if protein2cog[d] == "COG0085":
					id_hit2 = d +"|COG0085"
					range2 = protein2coords[id_hit2]
					r2 = range(range2[0], range2[1])
					meanloc2 = np.mean(range2)
								
					set1 = set(r1)
					inter = set1.intersection(r2)

					if int(len(inter)) > 50:
						protein2dups[id_hit2] = "false_hit"
					else:
						protein2dups[id_hit1] = "main_hit"
						protein2dups[id_hit2] = "secondary_hit"

						minrange = min(range1 + range2)
						maxrange = max(range1 + range2)
						protein2coords[id_hit1] = [minrange, maxrange]
									
						protein2align_length[id_hit1] = abs(maxrange - minrange)
						protein2length[rnap] = int(protein2length[rnap]) + int(protein2length[d])
						protein2score[rnap] = float(protein2score[rnap]) + float(protein2score[d])
						prot2protlist[rnap].append(d)
						prot2loc[rnap].append(meanloc2)
						num_proteins[id_hit1] +=1



			for item in protein2dups:
				#print item
				if protein2dups[item] == "main_hit" or protein2dups[item] == "single_besthit": #or protein2dups[item] == "NEXT" or protein2dups[item] == "SECO":
					items = item.split("|")
					protein = items[0]
					hit = items[1]

					protlist = prot2protlist[protein]
					loc_list = [str(loc) for loc in prot2loc[protein]]
					index_list = [i[0] for i in sorted(enumerate(loc_list), key=lambda x:x[1])]
					sorted_loc_list = [i[1] for i in sorted(enumerate(loc_list), key=lambda x:x[1])]

					sorted_prot_list = [protlist[index] for index in index_list]
					prot_str = ";".join(sorted_prot_list)

					loc_str = ";".join(sorted_loc_list)
					acc = protein2acc[protein]

					merged.write(protein +"\t"+ acc +"\t"+ taxon_dict[acc] +"\t"+ hit +"\t"+ protein2dups[item] +"\t"+ str(protein2length[protein]) +"\t"+ str(protein2score[protein]) +"\t"+ str(protein2align_length[item]) +"\t"+ str(num_proteins[item]) +"\t"+ prot_str +"\t"+ loc_str +"\n")
					#print sorted_prot_list

					if len(sorted_prot_list) > 1:
						tally = tally + len(sorted_prot_list)

						newrecord = SeqRecord(Seq("", IUPAC.protein), id=protein+" JOINED", name=protein+" JOINED", description=protein2acc[protein])
						for fragment in sorted_prot_list:
							#print fragment
							subrecord = seq_dict[fragment]
							subseq = subrecord.seq
							subseq = re.sub("\*", "", str(subseq))
							#print subseq
							#print record.seq
							#print type(subseq)
							newrecord.seq = newrecord.seq +""+ subseq

						#if len(newrecord.seq) > 800 and float(protein2score[protein]) > 400:
						final_proteins.append(newrecord)

					else:
						tally +=1
						#if len(seq_dict[protein].seq) > 800 and float(protein2score[protein]) > 400:
						final_proteins.append(seq_dict[protein])


SeqIO.write(final_proteins, merged_proteins, "fasta")
print tally




































