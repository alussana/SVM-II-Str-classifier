#!/usr/bin/env python3

import numpy as np
from sklearn import svm
import pickle
import gzip

def readDssp(dssp_dir, prot_id):
	dsspFile = open(dssp_dir + '/' + prot_id + '.dssp')
	dsspFile.readline()
	sequence = dsspFile.readline().rstrip()
	return(sequence)

def readPssm(pssmPath):
	dsspFile = open(pssmPath)
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
	linear_x = np.array(linear_x).astype(float)
	linear_x = linear_x.reshape(1,-1)
	return(linear_x)

def impute_II_str(pssm_dir, idList, out_dir, modelFile, windowSize):
	aa = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']	
	machine = pickle.load(gzip.open(modelFile))
	with open(idList) as ids:
		for prot_id in ids:
			print('Processing', prot_id)
			prot_id = prot_id.rstrip()
			pssm, aaDict = readPssm(pssm_dir + '/' + prot_id + '.pssm')
			dssp = ""
			for dssp_pos in range(len(pssm)):
				x = create_x_input(pssm, dssp_pos, windowSize, aa, aaDict)
				y = machine.predict(x)[0]
				if y == 0:
					dssp = dssp + 'H'
				elif y == 1:
					dssp = dssp + 'E'
				elif y == 2:
					dssp = dssp + '-'
			print_dssp(dssp, out_dir, prot_id)

def print_dssp(dssp, out_dir, prot_id):
	outFile = out_dir + '/' + prot_id + '.dssp'
	destination = open(outFile, 'w')
	print('>' + prot_id, file = destination)
	print(dssp, file = destination)
	destination.close()

if __name__ == '__main__':
	import argparse as ap
	parser = ap.ArgumentParser(description='Impute II Structure with SVM model')
	parser.add_argument('--pssm_dir', dest='pssm_dir', help='pssm directory')
	parser.add_argument('--id_list', dest='id_list', help='list of file names without extension')
	parser.add_argument('--model', dest='model_file', help='output of svm_training')
	parser.add_argument('--out_dir', dest='out_dir', help='output directory for predicted dssp files')
	parser.add_argument('--windowSize', dest='windowSize', type=int, help='length of the profile window')
	args = parser.parse_args()
	pssm_dir = args.pssm_dir
	idList = args.id_list
	out_dir = args.out_dir
	modelFile = args.model_file
	windowSize = args.windowSize
	impute_II_str(pssm_dir, idList, out_dir, modelFile, windowSize)
