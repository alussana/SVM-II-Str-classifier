#!/opt/conda/bin/python

import sys
import gzip

def filter_length(sourceFile, min_length):
	with gzip.open(sourceFile, 'rt') as fastaFile:
		for line in fastaFile:
			line = line.rstrip('\n')
			if line.startswith('>'):
				end_of_seq = 1
				header = line
			else:
				if len(line) >= min_length or end_of_seq == 0:
					if end_of_seq == 1:
						print(header)
					end_of_seq = 0
					print(line)

if __name__ == '__main__':
	sourceFile = sys.argv[1]
	min_length = int(sys.argv[2])
	filter_length(sourceFile, min_length)
