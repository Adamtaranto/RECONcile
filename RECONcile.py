#!/usr/bin/env python
#python 2.7.5 requires biopython
#RECONcile.py
#Version 1. Adam Taranto, March 2017
#Contact, Adam Taranto, adam.taranto@anu.edu.au

###########################################################################################
# Take clustered element fragment coordinates from RECON and write sequences to fasta for #
# alignment and consensus calling.                                                        #
###########################################################################################

import csv
import sys
import os
import subprocess
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def dochecks(args):
	# Check for reference fasa file
	if args.inFasta is None:
		sys.exit('No input fasta provided. Exiting.')
	# Check for eles file
	if args.eles is None:
		sys.exit('No RECON "eles" file provided. Exiting.')
	# Make outDir if does not exist else set to current dir
	if args.outDir:
		tempPathCheck(args)
		outDir = args.outDir
	else:
		outDir = os.getcwd() 
	return outDir

def tempPathCheck(args):
	absOutDir = os.path.abspath(args.outDir)
	if not os.path.isdir(absOutDir):
		os.makedirs(absOutDir)

def getLibrary(inFasta):
	#Populate dictionary with master set of seq records
	SeqMaster = dict()
	for seq_record in SeqIO.parse(inFasta, "fasta"):
		SeqMaster[seq_record.id] = str(seq_record.seq)
	return SeqMaster

def readEles(Cluster_file):
	with open(Cluster_file) as f:
		content = f.readlines()
	content = [x.strip().split() for x in content] 
	dictFam = dict()
	for row in content:
		fam = row[0]
		strand = row[2]
		name  = row[3]
		start = int(row[4]) - 1
		end   = int(row[5]) - 1
		if fam not in dictFam.keys():
			dictFam[fam] = list()
		dictFam[fam].append((strand,name,start,end))
	return dictFam

def chunkstring(string, length=80):
	return (string[0+i:length+i] for i in range(0, len(string), length))

def revComplement(seq):
    revcompl = lambda x: ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}[B] for B in x][::-1])
    return revcompl(seq)

def getFragments(clusterDict,famID,SeqMaster,errorHandle):
	for strand,name,start,end in clusterDict[famID]:
		try:
			SeqMaster[name]
		except:
			print('Reference seq not found: ' + name)
			errorHandle.write(name+"\n")
		else:
			if strand == "-1":
				FwdSeq = SeqMaster[name][start:end]
				seq = revComplement(FwdSeq)
				direction = "C"
			else:
				seq = SeqMaster[name][start:end]
				direction = "F"
			lable = "_".join([name,str(start),str(end),direction])
			yield (lable,seq)

def writeClusters(clusterDict, SeqMaster, outDir):
	# Open Error Log file
	errorPath = os.path.join(outDir, "errorlog.txt")
	errorHandle = open(errorPath,'w')
	# Write fragments for each cluster to output cluster file
	for famID in clusterDict.keys():
		outPath = os.path.join(outDir, "Cluster_" + famID + ".fa")
		outHandle = open(outPath,'w')
		for lable,seq in getFragments(clusterDict,famID,SeqMaster,errorHandle):
			fasta_name	= ">%s" % (lable)
			outHandle.write(fasta_name+"\n")
			for line in chunkstring(seq):
				outHandle.write(line+"\n")
		outHandle.close()
	errorHandle.close()

def main(args):
	# Check for required files + make output directory
	outDir = dochecks(args)
	# Make dictionary of Seq records keyed by sequence name
	SeqMaster = getLibrary(args.inFasta)
	# Read element names and fragment co-ords into dict keyed by family ID
	clusterDict = readEles(args.eles)
	#For each familes extract sequence fragments and write to file
	writeClusters(clusterDict, SeqMaster, outDir)


if __name__== '__main__':
	###Argument handling.
	parser = argparse.ArgumentParser(
		description='Extracts and orients TE frags belonging to RECON clusters.',
		prog='RECONcile')
	parser.add_argument("-i", "--inFasta",
		type=str,
		required=True,
		default= None,
		help="Multi fasta containing all TE sequences.")
	parser.add_argument("-e", "--eles",
		type=str,
		default= None,
		help="Space delimited ele file from RECON.")
	parser.add_argument("-d", "--outDir",
		type=str,
		default= None, 
		help="Directory for new cluster files to be written to.")

	args = parser.parse_args()

	main(args);