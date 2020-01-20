#!/usr/bin/env python3

import sys
import numpy as np

def statTools_read(input_file, separator = '\t'):
	data = {}
	for line in input_file:
		line = line.rstrip()
		values = line.split(separator)
		for i in range(len(values)):
			if i in data.keys():
				data[i].append(values[i])
			else:
				data[i] = [values[i]]
	table = [[] for i in range(len(data.keys()))]
	for i in range(len(table)):
		table[i] = data[i]
	data = np.array(table).astype(float)
	return(data)

def statTools_mean(data):
	m = list(map(np.mean, data))
	out = []
	for i in m:
		out.append(str(i))
	out = '\t'.join(out)
	print(out)

def statTools_stdDev(data):
	s = list(map(np.std, data))
	out = []
	for i in s:
		out.append(str(i))
	out = '\t'.join(out)
	print(out)

def statTools_stdErr(data):
	n = []
	out = []
	for i in data:
		n.append(len(data))
	n = np.array(n).astype(float)
	e = list(map(np.std, data))
	e = e / np.sqrt(n)
	for i in e:
		out.append(str(i))
	out = '\t'.join(out)
	print(out)

if __name__ == '__main__':
	import argparse as ap
	parser = ap.ArgumentParser(description='Python3 toolkit for statistics')
	parser.add_argument('--sep', help='Column separator')
	parser.add_argument('--mean', action='store_true', help='Compute the mean of each input column')
	parser.add_argument('--std_dev', action='store_true', help='Compute the standard deviation for each input column')
	parser.add_argument('--std_err', action='store_true', help='Compute the standard error for each input column')
	args = parser.parse_args()

	input_file = sys.stdin
	
	if args.sep:
		separator = args.sep
		data = statTools_read(input_file, separator)
	else:
		data = statTools_read(input_file)

	if args.mean:
		statTools_mean(data)
	if args.std_dev:
		statTools_stdDev(data)
	if args.std_err:
		statTools_stdErr(data)
