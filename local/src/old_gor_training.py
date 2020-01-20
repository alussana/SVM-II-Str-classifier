#!/opt/conda/bin/python

import numpy as np
import gzip

def readTable(fileName):
	#domainID = []
	sequence = []
	with gzip.open(fileName, 'rt') as sourceTable:
		for sourceLine in sourceTable:
			sourceInfo = sourceLine.split("\t")
			#domainID.append(sourceInfo[0])
			sequence.append(sourceInfo[1])
	return(sequence)

def computePropensity(windowSize, aa, fasta, dssp):
	H = [0] * windowSize
	E = [0] * windowSize
	c = [0] * windowSize
	for dssp_pos in range(len(dssp)):
		for x in range(len(fasta)):
		for win_pos in range(-windowSize // 2 + 1, windowSize // 2 + 1):	
			if dssp_pos + win_pos >= 0 and dssp_pos + win_pos < len(dssp):
				if dssp[x][i + windowSize // 2] == 'H':
				for j in range(windowSize):
					if fasta[x][i+j] == aa:
						H[j] = H[j] + 1
			elif dssp[x][i + windowSize // 2] == 'E':
				for j in range(windowSize):
					if fasta[x][i+j] == aa:
						E[j] = E[j] + 1
			i = i + 1
	return((H,E))

def generateTables(multitable):
	helixTable = []
	strandTable = []
	for aa in multitable:
		helixTable.append(aa[0])
		strandTable.append(aa[1])
	return(np.array(helixTable).astype(float), np.array(strandTable).astype(float))

def gor_training(fastaSource, dsspSource, outFile, windowSize):
	aa = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']
	fasta = readTable(fastaSource)
	dssp = readTable(dsspSource)
	propensity17 = []
	for aminoAcid in aa:
		propensity17.append(computePropensity(windowSize, aminoAcid, fasta, dssp))
	helixTable, strandTable = generateTables(propensity17)
	tot_counts = 0
	for j in range(len(aa)):
		tot_counts = tot_counts + helixTable[j][0]
	for i in range(windowSize):
		for j in range(len(aa)):
			helixTable[j][i] = helixTable[j][i] / tot_counts
			strandTable[j][i] = strandTable[j][i] / tot_counts
	return(helixTable, strandTable)

def print_training_file(helixTable, strandTable, aa, windowSize, outFile):
	output = open(outFile, 'w')
	positions = [str(i) for i in range(windowSize)]
	first_line = '#\t' + '\t'.join(positions)
	print('>H', file = output)
	print(first_line, file = output)
	for j in range(len(aa)):
		line = aa[j]
		for i in range(len(positions)):
			line = line + '\t%f'%(helixTable[j][i])
		print(line, file = output)
	print('>E', file = output)
	print(first_line, file = output)
	for j in range(len(aa)):
		line = aa[j]
		for i in range(len(positions)):
			line = line + '\t%f'%(strandTable[j][i])
		print(line, file = output)

if __name__ == '__main__':
	import argparse as ap
	parser = ap.ArgumentParser(description='Train GOR method with input data.')
	parser.add_argument('--fasta', dest='fasta_source', help='fasta file')
	parser.add_argument('--dssp', dest='dssp_source', help='dssp file')
	parser.add_argument('--out', dest='out_file', help='output file name')
	parser.add_argument('--window_size', dest='window_size', type=int, help='max number of residues that define the training windows size')
	args = parser.parse_args()
	aa = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']
	# fastaSource and dsspSource have to be consistently sorted
	fastaSource = args.fasta_source
	dsspSource = args.dssp_source
	outFile = args.out_file
	windowSize = args.window_size
	helixTable, strandTable = gor_training(fastaSource, dsspSource, outFile, windowSize)
	print_training_file(helixTable, strandTable, aa, windowSize, outFile)
