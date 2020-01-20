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

def find_all_zero_pssm(pssm_dir, idList):
	all_zero_ids = []
	with open(idList) as ids:
		for prot_id in ids:
			prot_id = prot_id.rstrip()
			pssm, aaDict = readPssm(pssm_dir, prot_id)
			## check for all-zero matrices, discard them, print is STDOUT the prot_id:
			check_pos = []
			for i in range(len(pssm)):
				check_aa = pssm[i]
				c = check_aa == 0
				it = iter(c)
				check_pos.append(all(it))
			it = iter(c)
			if all(it):
				all_zero_ids.append(prot_id)		
	return(all_zero_ids)

def print_list(l):
	for i in l:
		i = i.rstrip()
		print(i)

if __name__ == '__main__':
	import argparse as ap
	parser = ap.ArgumentParser(description='Reports pssm ids for those having zero vectors in all positions')
	parser.add_argument('--pssm_dir', dest='pssm_dir', help='pssm directory')
	#parser.add_argument('--dssp_dir', dest='dssp_dir', help='dssp directory')
	parser.add_argument('--id_list', dest='id_list', help='list of file names without extension')
	#parser.add_argument('--out', dest='out_file', help='output file name')
	args = parser.parse_args()
	pssm_dir = args.pssm_dir
	#dssp_dir = args.dssp_dir
	idList = args.id_list
	#outFile = args.out_file

	all_zero_ids = find_all_zero_pssm(pssm_dir, idList)
	print_list(all_zero_ids)
