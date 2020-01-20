#!/opt/conda/bin/python

import numpy as np

def readPssm(pssm_dir, prot_id):
	dsspFile = open(pssm_dir + '/' + prot_id + '.pssm')
	dsspInfo = []
	for line in dsspFile:
		dsspInfo.append(line.strip())
	dsspFile.close()
	aa = dsspInfo[2].split()[:20]
	aaDict = dict(zip(aa, list(range(20))))
	pssmValues = []
	for position in range(3, len(dsspInfo) - 6):
		dsspLine = dsspInfo[position].split()[2:]
		dsspLine = dsspLine[20:40]
		dsspLine = [float(i) / 100 for i in dsspLine]
		pssmValues.append(dsspLine)
	pssmValues = np.array(pssmValues).astype(float)
	return(pssmValues, aaDict)

def readDssp(dssp_dir, prot_id):
	dsspFile = open(dssp_dir + '/' + prot_id + '.dssp')
	dsspFile.readline()
	sequence = dsspFile.readline().rstrip()
	return(sequence)

def update_propensieties(table, priorsTable, pssmValues, dssp_pos, windowSize, aa, aaDict):
	for win_pos in range(-windowSize // 2 + 1, windowSize // 2 + 1):
		if dssp_pos + win_pos >= 0 and dssp_pos + win_pos < len(pssmValues):
			propensities = pssmValues[dssp_pos + win_pos]
			for i in range(len(aa)):
				aminoAcid = aa[i]
				aa_index = aaDict[aminoAcid]
				table[i][win_pos + windowSize // 2] = table[i][win_pos + windowSize // 2] + propensities[aa_index]
				priorsTable[i][win_pos + windowSize // 2] = priorsTable[i][win_pos + windowSize // 2] + propensities[aa_index]
	return(table, priorsTable)

def gor_train(pssm_dir, dssp_dir, outFile, windowSize, idList):
	aa = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']
	helixTable = np.array([[0] * windowSize] * 20).astype(float)
	strandTable = np.array([[0] * windowSize] * 20).astype(float)
	coilTable = np.array([[0] * windowSize] * 20).astype(float)
	priorsTable = np.array([[0] * windowSize] * 20).astype(float)
	with open(idList) as ids:
		n_profiles = 0
		n_residues = 0
		for prot_id in ids:
			prot_id = prot_id.rstrip()
			pssmValues, aaDict = readPssm(pssm_dir, prot_id)
			dssp = readDssp(dssp_dir, prot_id)
			if len(pssmValues) != len(dssp):
				print("E: %s has length inconsistency"%(prot_id))
			else:
				n_profiles = n_profiles + 1
				n_residues = n_residues + len(dssp)
				for dssp_pos in range(len(dssp)):
					if dssp[dssp_pos] == 'H':
						helixTable, priorsTable = update_propensieties(helixTable, priorsTable, pssmValues, dssp_pos, windowSize, aa, aaDict)	
					elif dssp[dssp_pos] == 'E':
						strandTable, priorsTable = update_propensieties(strandTable, priorsTable, pssmValues, dssp_pos, windowSize, aa, aaDict)
					elif dssp[dssp_pos] == '-':					
						coilTable, priorsTable = update_propensieties(coilTable, priorsTable, pssmValues, dssp_pos, windowSize, aa, aaDict)
					else:
						print("E: %s: unexpected character in dssp file"%(prot_id))
				#helixTable = helixTable / len(dssp)
				#strandTable = strandTable / len(dssp)
				#coilTable = coilTable / len(dssp)
				#priorsTable = priorsTable / len(dssp)
	helixTable = helixTable / n_residues
	strandTable = strandTable / n_residues
	coilTable = coilTable / n_residues
	priorsTable = priorsTable / n_residues
	## compute information values
	H_freq = 0
	for i in range(len(helixTable)):
		H_freq = H_freq + helixTable[i][windowSize // 2]
	E_freq = 0
	for i in range(len(strandTable)):
		E_freq = E_freq + strandTable[i][windowSize // 2]
	c_freq = 0
	for i in range(len(coilTable)):
		c_freq = c_freq + coilTable[i][windowSize // 2]
	for amino in range(len(helixTable)):
		for pos in range(len(helixTable[0])):
			helixTable[amino][pos] = np.log2(helixTable[amino][pos] / (priorsTable[amino][pos] * H_freq))
			strandTable[amino][pos] = np.log2(strandTable[amino][pos] / (priorsTable[amino][pos] * E_freq))
			coilTable[amino][pos] = np.log2(coilTable[amino][pos] / (priorsTable[amino][pos] * c_freq))
	print('Training was performed on %i protein profiles'%n_profiles)
	return(helixTable, strandTable, coilTable, priorsTable)

def print_training_file(helixTable, strandTable, coilTable, priorsTable, aa, windowSize, outFile):
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
	print('>c', file = output)
	print(first_line, file = output)
	for j in range(len(aa)):
		line = aa[j]
		for i in range(len(positions)):
			line = line + '\t%f'%(coilTable[j][i])
		print(line, file = output)
	print('>priors', file = output)
	for j in range(len(aa)):
		line = aa[j]
		for i in range(len(positions)):
			line = line + '\t%f'%(priorsTable[j][i])
		print(line, file = output)

if __name__ == '__main__':
	## parse arguments
	import argparse as ap
	parser = ap.ArgumentParser(description='Train GOR method with input data.')
	parser.add_argument('--pssm_dir', dest='pssm_dir', help='pssm directory')
	parser.add_argument('--dssp_dir', dest='dssp_dir', help='dssp directory')
	parser.add_argument('--id_list', dest='id_list', help='list of file names without extension')
	parser.add_argument('--out', dest='out_file', help='output file name')
	parser.add_argument('--window_size', dest='window_size', type=int, help='max number of residues that define the training windows size')
	args = parser.parse_args()
	pssm_dir = args.pssm_dir
	dssp_dir = args.dssp_dir
	idList = args.id_list
	outFile = args.out_file
	windowSize = args.window_size
	helixTable, strandTable, coilTable, priorsTable = gor_train(pssm_dir, dssp_dir, outFile, windowSize, idList)
	aa = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']
	print_training_file(helixTable, strandTable, coilTable, priorsTable, aa, windowSize, outFile)
