#!/opt/conda/bin/python

import sys
import gzip

def get_dssp(sourceName, fastaName, dsspName):
	fastaFile = open(fastaName, 'a')
	dsspFile = open(dsspName, 'a')
	with gzip.open(sourceName, 'rt') as sourceFile:
		seq = ""
		for line in sourceFile:
			if line.startswith('>'):
				if seq != "":
					print(seq, file = fastaFile)
					print(dssp, file = dsspFile)
				chain_id = line.rstrip()
				chain_ended = 1
				seq = ""
				dssp = ""
				print(chain_id, file = fastaFile)
				print(chain_id, file = dsspFile)
			elif line.startswith('  #'):
				chain_ended = 0
			elif chain_ended == 0:
				aa = line[13]
				structure = line[16]
				if structure == 'G' or structure == 'I':
					structure = 'H'
				elif structure == 'B':
					structure = 'E'
				elif structure == 'T' or structure == 'S' or structure == ' ':
					structure = '-'
				seq = seq + aa
				dssp = dssp + structure
	fastaFile.close()
	dsspFile.close()

if __name__ == '__main__':
	sourceName = sys.argv[1]
	fastaName = sys.argv[2]
	dsspName = sys.argv[3]
	get_dssp(sourceName, fastaName, dsspName)
