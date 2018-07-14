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

marker_out = "/home/frankaylward/Desktop/marker_gene_benchmarking/genbank_curated_set/hmm_out/speci/"
downloads = "/home/frankaylward/Desktop/marker_gene_benchmarking/genbank_curated_set/genomes/"
info = open("/home/frankaylward/Desktop/marker_gene_benchmarking/genbank_curated_set/genome_manifest_partial_2.txt", "r")

merged = open("merged_COG0085_raw.txt", "w")
merged.write("protein\tacc\thit\tspecies\tdomain\tphylum\tcategory\tlength\tscore\talign_length\tnum_prots\tprot_ids\tprotein_locs\n")

merged_proteins = open("merged_COG0085.faa", "w")
final_proteins = []

prox = int(10)

# genome list
genomes = open("good_quality_genomes.list", "r")
genome_dict = {}
for j in genomes.readlines():
	line = j.rstrip()
	genome_dict[line] = line


# get taxonomic info
taxon_dict = {}
for j in info:
	line = j.rstrip("\n")
	tabs = line.split("\t")
	name = tabs[0]
	taxonomic = "\t".join(tabs[2:5])
	taxon_dict[name] = taxonomic


marker_tally = defaultdict(int)
exceptions = open("exceptions.txt", "w")
tally = 0

for i in os.listdir(downloads):
	if i in genome_dict:
		
		###############################################
		########## Part 1: Parse the GFF file #########
		###############################################
		region2proteins = defaultdict(list)
		protein2region = {}

		# get gff dictionaries for this accession
		acc = re.sub(".speci_out.parsed", "", i)
		path = os.path.join(downloads, i)
		print i

		path_to_gff = glob.glob(path+"/*.prodigal.gff")
		if len(path_to_gff) < 1:
			path_to_gff = glob.glob(path+"/*genomic.gff")
			#print path_to_gff
			gff_handle = open(path_to_gff[0], "r")

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
						for m in semi:
							if "protein_id" in m:
								prot = re.sub("protein_id=", "", m)

						region2proteins[region].append(prot)
						protein2region[prot] = region

		else:
			#print path, path_to_gff
			gff_handle = open(path_to_gff[0], "r")

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

		path_to_faa = glob.glob(path+"/*.faa")
		seq_dict = SeqIO.to_dict(SeqIO.parse(path_to_faa[0], "fasta"))

		################################################################
		########## Part 3: Get alignment coords from domtblout #########
		################################################################
		marker_file = os.path.join(marker_out, i +".speci_dom_out.parsed")
		parsed = open(os.path.join(marker_out, marker_file), "r")
		done = {}
		protein2coords = defaultdict(list)
		protein2align_length = {}

		for n in parsed.readlines():
			line = n.rstrip()
			tabs = line.split("\t")
			protein = tabs[0]
			annot = tabs[1]
			start = int(tabs[3])
			end = int(tabs[4])
			id_hit = protein +"|"+ annot

			#print protein, annot, start, end, id_hit
			protein2coords[id_hit].append(start)
			protein2coords[id_hit].append(end)

			align_length = abs(end - start)
			protein2align_length[id_hit] = align_length


		#########################################################
		########## Part 4: Parse the marker Annote file #########
		#########################################################
		marker_file = os.path.join(marker_out, i +".speci_out.parsed")
		parsed = open(os.path.join(marker_out, marker_file), "r")
		done = {}
		protein2cog = {}
		protein2acc = {}
		protein2score = {}
		protein2category = {}
		protein2length = {}
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
				acc = re.sub(".speci_out.parsed", "", i)
				nr = acc +"_"+ annot

				acc2hits[acc]

				protein2cog[protein] = annot
				protein2acc[protein] = acc
				protein2score[protein] = hmm_score
				protein2length[protein] = prot_length

				if nr in done:
					category = "NH"
				else:
					category = "BH"
					done[nr] = annot
					protein2dups[id_hit]

				protein2category[protein] = category
		parsed.close()

		parsed = open(os.path.join(marker_out, marker_file), "r")
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

					if index + prox > length:
						right_bound = length
					else:
						right_bound = index + prox

					if index - prox < 0:
						left_bound = 0
					else:
						left_bound = index - prox

					orf_set = set(protein_list[left_bound:right_bound])

					#index = int(protein_list.index(protein))
					#length = len(protein_list)

					#orf_set = set(protein_list)

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

								#if protein == "KZX16681.1":
								#	print range1, range2, inter, region

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
			if protein2dups[item] == "MAIN" or protein2dups[item] == "PRIMARY" or protein2dups[item] == "NEXT" or protein2dups[item] == "SECO":
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














