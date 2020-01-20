#!/usr/bin/env python3

import numpy as np
from sklearn import svm

def readDssp(dssp_dir, prot_id):
	dsspFile = open(dssp_dir + '/' + prot_id + '.dssp')
	dsspFile.readline()
	sequence = dsspFile.readline().rstrip()
	return(sequence)

def readPssm(pssm_dir, prot_id):
	dsspFile = open(pssm_dir + '/' + prot_id + '.pssm')
	dsspInfo = []
	for line in dsspFile:
		dsspInfo.append(line.strip())
	dsspFile.close()
	aa = dsspInfo[2].split()[:20]
	aaDict = dict(zip(aa, list(range(20))))
	pssm = []
	for position in range(3, len(dsspInfo) - 6):
		dsspLine = dsspInfo[position].split()[2:]
		dsspLine = dsspLine[20:40]
		dsspLine = [float(i) / 100 for i in dsspLine]
		pssm.append(dsspLine)
	pssm = np.array(pssm).astype(float)
	return(pssm, aaDict)

def convert_dssp(dssp, prot_id):
	y = []
	for i in range(len(dssp)):
		if dssp[i] == 'H':
			y.append(0)
		elif dssp[i] == 'E':
			y.append(1)
		elif dssp[i] == '-':
			y.append(2)
		else:
			print("E: %s: dssp contains unexpected characters"%(prot_id))
	return(y)

def create_x_input(pssm, dssp_pos, windowSize, aa, aaDict):
	x = [([0] * windowSize) for i in range(len(aa))]
	for win_pos in range(-windowSize // 2 + 1, windowSize // 2 + 1):
		if dssp_pos + win_pos >= 0 and dssp_pos + win_pos < len(pssm):
			propensities = pssm[dssp_pos + win_pos]
			for i in range(len(aa)):
				aminoAcid = aa[i]
				aa_index = aaDict[aminoAcid]
				x[i][win_pos + windowSize // 2] = x[i][win_pos + windowSize // 2] + propensities[aa_index]
	linear_x = []
	for aa_row in x:
		for value in aa_row:
			linear_x.append(value)
	return(linear_x)

def build_training_examples(pssm_dir, dssp_dir, idList, windowSize):
	aa = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']
	x = []
	y = []
	with open(idList) as ids:
		n_profiles = 0
		n_residues = 0
		for prot_id in ids:
			prot_id = prot_id.rstrip()
			pssm, aaDict = readPssm(pssm_dir, prot_id)
			dssp = readDssp(dssp_dir, prot_id)
			if len(pssm) != len(dssp):
				print("E: %s has length inconsistency"%(prot_id))
			else:
				n_profiles = n_profiles + 1
				n_residues = n_residues + len(dssp)
				dssp = convert_dssp(dssp, prot_id)
				## check for all-zero matrices, discard them, print is STDOUT the prot_id:
				check_pos = []
				for i in range(len(pssm)):
					check_aa = pssm[i]
					c = check_aa == 0
					it = iter(c)
					check_pos.append(all(it))
				it = iter(c)
				if all(it):
					print('%s\tall-zero pssm is discarded from training'%(prot_id))
				else:
					for dssp_pos in range(len(dssp)):
						new_x_input = create_x_input(pssm, dssp_pos, windowSize, aa, aaDict)
						x.append(new_x_input)
						y.append(dssp[dssp_pos])
	x = np.array(x).astype(float)
	y = np.array(y).astype(int)
	#print("Training data was retrieved on %i protein profiles including $i residues"%(n_profiles, n_residues))
	return(x, y)

def SVM_train(pssm_dir, dssp_dir, idList, windowSize, outFile, kernel, C, gamma):
	x, y = build_training_examples(pssm_dir, dssp_dir, idList, windowSize)
	machine = svm.SVC(C = C, kernel = kernel, gamma = gamma)
	machine.fit(x, y)
	return(machine)

def export_SVM(machine, outFile):
	import pickle, gzip
	pickle.dump(machine, gzip.open(outFile, 'w'))

if __name__ == '__main__':
	import argparse as ap
	parser = ap.ArgumentParser(description='Train a SVM to impute secondary structure of proteins')
	parser.add_argument('--pssm_dir', dest='pssm_dir', help='pssm directory')
	parser.add_argument('--dssp_dir', dest='dssp_dir', help='dssp directory')
	parser.add_argument('--id_list', dest='id_list', help='list of file names without extension')
	parser.add_argument('--out', dest='out_file', help='output file name')
	parser.add_argument('--window_size', dest='window_size', type=int, help='max number of residues that define the training windows size')
	parser.add_argument('--kernel', dest='kernel', type=str, help='kernel type for SVM')
	parser.add_argument('--C', dest='C', type=float, help='C parameter for SVM')
	parser.add_argument('--gamma', dest='gamma', type=float, help='gamma parameter for SVM')
	args = parser.parse_args()
	pssm_dir = args.pssm_dir
	dssp_dir = args.dssp_dir
	idList = args.id_list
	outFile = args.out_file
	windowSize = args.window_size
	kernel = args.kernel
	C = args.C
	gamma = args.gamma
	aa = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']
	machine = SVM_train(pssm_dir, dssp_dir, idList, windowSize, outFile, kernel, C, gamma)
	export_SVM(machine, outFile)
