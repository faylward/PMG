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

#marker_out = "/home/frankaylward/Desktop/marker_gene_benchmarking/genbank_curated_set/hmm_out/speci/"
downloads = "/home/frankaylward/Desktop/marker_gene_benchmarking/metagenomes/uba_assemblies"

merged = open("merged_COG0085.txt", "w")
merged.write("protein\tacc\thit\tcategory\tlength\tscore\talign_length\tnum_prots\tprot_ids\tprotein_locs\n")

merged_proteins = open("merged_COG0085.faa", "w")
final_proteins = []

#genomes = open("good_quality_genomes.list", "r")
#genome_dict = {}
#for j in genomes.readlines():
#	line = j.rstrip()
#	genome_dict[line] = line


marker_tally = defaultdict(int)
exceptions = open("exceptions.txt", "w")
tally = 0
placeholder = 1

path_to_gff = glob.glob(downloads+"/*.gff")
for i in path_to_gff:

	###############################################
	######### Part 1: Parse the GFF file ##########
	###############################################
	region2proteins = defaultdict(list)
	protein2region = {}

	# get gff dictionaries for this accession
	acc1 = re.sub("_contigs.gff", "", i)
	acc = re.sub(downloads, "", acc1)
	acco = re.sub("/", "", acc)

	#print path, path_to_gff
	print acco
	gff_handle = open(i, "r")

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

	##################################################################
	########## Part 2: Get protein sequences in a dictionary #########
	##################################################################

	#path_to_faa = glob.glob(path+"/*.faa")
	path_to_faa = re.sub(".gff", ".faa", i)
	seq_dict = SeqIO.to_dict(SeqIO.parse(path_to_faa, "fasta"))

	################################################################
	########## Part 3: Get alignment coords from domtblout #########
	################################################################
	#marker_file = os.path.join(marker_out, i +".speci_dom_out.parsed")
	marker_file = re.sub(".gff", ".hmm_dom_out.parsed", i)
	parsed = open(marker_file, "r")
	done = {}
	protein2coords = defaultdict(list)
	protein2length = {}
	protein2align_length = {}

	for n in parsed.readlines():
		line = n.rstrip()
		tabs = line.split("\t")
		protein = tabs[0]
		annot = tabs[1]
		start = int(tabs[4])
		end = int(tabs[5])
		id_hit = protein +"|"+ annot
		protein2length[protein] = tabs[2]

		#print protein, annot, start, end, id_hit
		protein2coords[id_hit].append(start)
		protein2coords[id_hit].append(end)

		align_length = abs(end - start)
		protein2align_length[id_hit] = align_length

	#########################################################
	########## Part 4: Parse the marker Annote file #########
	#########################################################
	#marker_file = os.path.join(marker_out, i +".hmm_out.parsed")
	marker_file = re.sub(".gff", ".hmm_out.parsed", i)
	parsed = open(marker_file, "r")
	done = {}
	protein2cog = {}
	protein2acc = {}
	protein2score = {}
	protein2category = {}

	protein2dups = defaultdict(lambda:"PRIMARY")
	acc2hits = defaultdict(lambda:"MAIN")

	#print i
	for n in parsed.readlines():
		line = n.rstrip()
		tabs = line.split("\t")
		protein = tabs[0]
		annot = tabs[2]

		if annot == "COG0085":

			id_hit = protein +"|"+ annot
			hmm_score = tabs[3]
			prot_length = tabs[4]
			acc = re.sub(".hmm_out.parsed", "", i)
			nr = acc +"_"+ annot

			acc2hits[acc]

			protein2cog[protein] = annot
			protein2acc[protein] = acco
			protein2score[protein] = hmm_score

			if nr in done:
				category = "NH"
			else:
				category = "BH"
				done[nr] = annot
				protein2dups[id_hit]

			protein2category[protein] = category
	parsed.close()

	parsed = open(marker_file, "r")
	#for protein in protein2cog:
	already_done = []
	num_proteins = defaultdict(lambda:int(1))
	prot2protlist = defaultdict(list)
	prot2loc = defaultdict(list)
		
	for n in parsed.readlines():
		line = n.rstrip()
		tabs = line.split("\t")
		protein = tabs[0]
		annot = tabs[2]
		#print protein, annot

		if annot == "COG0085":

			if annot in already_done:
				if protein2dups[protein +"|"+ annot] == "SECO":
					pass
				else:
					protein2dups[protein +"|"+ annot] = "NEXT"
			else:
				already_done.append(annot)
				prot2protlist[protein].append(protein)

				id_hit1 = protein +"|"+ annot
				range1 = protein2coords[id_hit1]
				r1 = range(range1[0], range1[1])
				meanloc1 = np.mean(range1)
				prot2loc[protein].append(meanloc1)

				region = protein2region[protein]

				#protein2dups[protein +"|"+ annot] = "NEXT"

				protein_list = region2proteins[region]

				index = int(protein_list.index(protein))
				length = len(protein_list)

				orf_set = set(protein_list)

				if protein in orf_set:
					orf_set.remove(protein)
				else:
					exceptions.write(start +"\t"+ end +"\t"+ left_bound +"\t"+ right_bound +"\t"+ region_length +"\n")

				for d in orf_set:
					if d in protein2cog:
						if protein2cog[protein] == protein2cog[d]:
							cog = protein2cog[protein]
							#id_hit1 = protein +"|"+ cog
							id_hit2 = d +"|"+ cog

							range2 = protein2coords[id_hit2]
							r2 = range(range2[0], range2[1])
							meanloc2 = np.mean(range2)
								
							set1 = set(r1)
							inter = set1.intersection(r2)

							if int(len(inter)) > 10:
								protein2dups[id_hit2] = "NEXT"
								pass
							else:
								protein2dups[id_hit1] = "MAIN"
								protein2dups[id_hit2] = "SECO"

								minrange = min(range1 + range2)
								maxrange = max(range1 + range2)
								protein2coords[id_hit1] = [minrange, maxrange]
									
								protein2align_length[id_hit1] = abs(maxrange - minrange)
								protein2length[protein] = int(protein2length[protein]) + int(protein2length[d])
								protein2score[protein] = float(protein2score[protein]) + float(protein2score[d])
								prot2protlist[protein].append(d)
								prot2loc[protein].append(meanloc2)
								num_proteins[id_hit1] +=1

	for item in protein2dups:
		if protein2dups[item] == "MAIN" or protein2dups[item] == "PRIMARY" or protein2dups[item] == "NEXT":
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


			merged.write(acco +"_"+ protein +"\t"+ protein2acc[protein] +"\t"+ hit +"\t"+ protein2dups[item] +"\t"+ str(protein2length[protein]) +"\t"+ str(protein2score[protein]) +"\t"+ str(protein2align_length[item]) +"\t"+ str(num_proteins[item]) +"\t"+ prot_str +"\t"+ loc_str +"\n")
			#print sorted_prot_list

			if len(sorted_prot_list) > 1:
				tally = tally + len(sorted_prot_list)

				newrecord = SeqRecord(Seq("", IUPAC.protein), id=acco+"_"+protein+"_JOINED", name=acco+"_"+protein+"_JOINED", description=protein2acc[protein])
				for fragment in sorted_prot_list:
					#print fragment
					subrecord = seq_dict[fragment]
					subseq = subrecord.seq
					subseq = re.sub("\*", "", str(subseq))
					#print subseq
					#print record.seq
					#print type(subseq)
					newrecord.seq = newrecord.seq +""+ subseq

				#if len(newrecord.seq) > 900 :# and protein2score[protein] > 500:
				final_proteins.append(newrecord)

			else:
				tally +=1
				#if len(seq_dict[protein].seq) > 900: #and protein2score[protein] > 500:
				record = seq_dict[protein]
				record.id = acco+"_"+record.id
				record.name = acco+"_"+record.id
				final_proteins.append(record)


SeqIO.write(final_proteins, merged_proteins, "fasta")
print tally














