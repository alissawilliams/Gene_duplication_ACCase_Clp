# local_blast5.py
# script to loop through results of a local blast search 
# and use them in reciprocal blasts
# any sequences that pass the reciprocal blast are put into a final fasta file
# Alissa Williams
# June 16, 2020

# to run, include:
# 1) name of database (subject of blast), with path and extension
# 2) name of the final file, with path and extension
# 3) name of the xml blast results file with path and extension


from Bio import SeqIO
from Bio.Blast import NCBIXML
import os
import sys

# use parameters from script call
databasename = sys.argv[1]
finalfilename = sys.argv[2]
xmlfilename = sys.argv[3]

# read in fasta database as a dictionary
seqdict = SeqIO.to_dict(SeqIO.parse(databasename, "fasta"))

# file for hits that are actually clpR2
finalfile = open(finalfilename, "w")

# open up xml output file
blast_records = NCBIXML.parse(open(xmlfilename))
# loop through blast results
for blast_record in blast_records:
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			# get name of sequence
			fullTitle = alignment.title
			# split name into actual name found in fasta file
			actualName = fullTitle.split(" ",2)[1]
			# find the sequence with that name 
			fullseq = seqdict.get(actualName).seq
			# write the output to a file
			newfile = open("currenthit_clpR2.fas", "w")
			newfile.write(">" + actualName + "\n")
			newfile.write(str(fullseq) + "\n")
			newfile.close()
			# reciprocal blast
			#make command
			recipcommand = "blastp -query currenthit_clpR2.fas -db Athaliana_167_TAIR10.protein.fa -out currenthit_recipblast_clpR2.txt"
			# actually run command, output will be a text file
			os.system(recipcommand)
			# get name of first hit using grep 
			twolinesabove = "Sequences producing significant alignments:                          (Bits)  Value"
			firsthitcommand = "grep -A 2 " + "\"" + twolinesabove + "\"" + " currenthit_recipblast_clpR2.txt" + " | tail -n 1"
			tophit = os.popen(firsthitcommand).read()
			tophit = tophit.split(" ")[2]
			tophit = tophit.split(".")[0]
			#print(tophit)
			# write the Amborella sequence to the final file if the top reciprocal blast hit is clpR2
			if str(tophit) == "AT1G12410":
				#print("entering if statement\n")
				finalfile.write(">" + actualName + "\n")	
				finalfile.write(str(fullseq) + "\n")
				#print(hitseq)

finalfile.close()








