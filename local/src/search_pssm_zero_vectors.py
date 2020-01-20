#!/usr/bin/env python3

import numpy as np

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

def search_zero_vectors(pssm_dir, idList):
	zero_vectors_ids = []
	counts = []
	with open(idList) as ids:
		for prot_id in ids:
			count = 0
			prot_id = prot_id.rstrip()
			pssm, aaDict = readPssm(pssm_dir, prot_id)
			for pos in pssm:
				if all(pos == 0):
					count = count + 1
			#if count != 0:
			#	zero_vectors_ids.append(prot_id)
			#	counts.append(count)
			zero_vectors_ids.append(prot_id)
			counts.append(count)
	return(zero_vectors_ids, counts)

def print_list_2col(l0, l1):
	if len(l0) != len(l1):
		print('E: Lists have different length')
		return
	else:
		for i in range(len(l0)):
			print(l0[i], l1[i])

if __name__ == '__main__':
	import argparse as ap
	parser = ap.ArgumentParser(description='Reports pssm ids for thos having at least one zero vector, and the number of zero vectors')
	parser.add_argument('--pssm_dir', dest='pssm_dir', help='pssm directory')
	parser.add_argument('--id_list', dest='id_list', help='list of file names without extension')
	args = parser.parse_args()
	pssm_dir = args.pssm_dir
	idList = args.id_list
	
	# ids: presence of zero vectors in the pssm that passed check_pssm.py
	ids, counts = search_zero_vectors(pssm_dir, idList)
	print_list_2col(ids, counts)
