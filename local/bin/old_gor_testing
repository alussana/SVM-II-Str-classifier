#!/opt/conda/bin/python

import numpy as np

def read_model(modelFile):
	aa = []
	with open(modelFile) as modelSource:
		aa = []
		H = 0
		E = 0
		c = 0
		for line in modelSource:
			line = line.rstrip()
			if line.startswith('>'):
				if line[1] == 'H':
					H = 1;	E = 0; c = 0
					helixTable = []
				elif line[1] == 'E':
					H = 0; E = 1; c = 0
					strandTable = []
				elif line[1] == 'c':
					H = 0; E = 0; c = 1
					coilTable = []
			elif  H == 1 and E == 0 and c == 0 and line.startswith('#') == False:
				info = line.split('\t')
				aa.append(info[0])
				helixTable.append(np.array(info[1:]).astype(float))
			elif H == 0 and E == 1 and c == 0 and line.startswith('#') == False:
				info = line.split('\t')
				strandTable.append(np.array(info[1:]).astype(float))
			elif H == 0 and E == 0 and c == 1 and line.startswith('#') == False:
				coilTable.append(np.array(info[1:]).astype(float))
		windowSize = len(helixTable[0] - 1)
	return(helixTable, strandTable, coilTable, aa, windowSize)

def impute_II_str(helixTable, strandTable, coilTable, aa, windowSize, fastaFile, outFile):
	with open(fastaFile) as fastaSource:
		header = fastaSource.readline().rstrip()
		sequence = ''
		for line in fastaSource:
			sequence = sequence + line.rstrip()
		II_structure = ''
		for i in range(len(sequence)):
			H_propensity = 0
			E_propensity = 0
			c_propensity = 0
			for j in range(-windowSize // 2, windowSize // 2 + 1):
				if i - j >= 0:
					H_propensity = H_propensity + helixTable[aa_dict[sequence[i]]][j]
					E_propensity = E_propensity + strandTable[aa_dict[sequence[i]]][j]
			if H_propensity > c_propensity:
				if H_propensity > E_propensity:
					II_structure = II_structure + 'H'
				else:
					II_structure = II_structure + 'E'
			else:
				if c_propensity > E_propensity:
					II_structure = II_structure + 'c'
				else:
					II_structure = II_structure + 'E'
	output = open(outFile, 'w')
	print(header, file=output)
	print(sequence, file=output)
	close(output)

if __name__ == '__main__':
	import argparse as ap
	parser = ap.ArgumentParser(description='Train GOR method with input data.')
	parser.add_argument('--input', dest='fasta_source', help='fasta file')
	parser.add_argument('--gor_model', dest='model_file', help='output of gor_training')
	parser.add_argument('--out', dest='out_file', help='output file name')
	args = parser.parse_args()
	fastaFile = args.fasta_source
	modelFile = args.model_file
	outFile = args.out_file
	helixTable, strandTable, coilTable, aa, windowSize = read_model(modelFile)
	impute_II_str(helixTable, strandTable, coilTable, aa, windowSize, fastaFile, outFile)
