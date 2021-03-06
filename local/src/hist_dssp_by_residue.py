#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

def getCounts(fileName):
	aa = []
	H = []
	E = []
	c = []
	with open(fileName) as sourceTable:
		for sourceLine in sourceTable:
			sourceInfo = sourceLine.rstrip().split("\t")
			aa.append(sourceInfo[0])
			H.append(int(sourceInfo[1]))
			E.append(int(sourceInfo[2]))
			c.append(int(sourceInfo[3]))
	return(aa, H, E, c)

def getFrequencies(aa, H, E, c):
	tot_H = 0
	tot_E = 0
	tot_c = 0
	tot = []
	for i in range(len(aa)):
		tot_H = tot_H + H[i]
		tot_E = tot_E + E[i]
		tot_c = tot_c + c[i]
		tot.append(H[i] + E[i] + c[i])
	tot_overall = tot_H + tot_E + tot_c
	freq_H = []
	freq_E = []
	freq_c = []
	freq_tot = []
	for i in range(len(aa)):
		freq_H.append(H[i]/tot_H)
		freq_E.append(E[i]/tot_E)
		freq_c.append(c[i]/tot_c)
		freq_tot.append(tot[i]/tot_overall)
	return(aa, freq_tot, freq_H, freq_E, freq_c)

## TODO https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/barchart.html#sphx-glr-gallery-lines-bars-and-markers-barchart-py
def multybarsHist(outFile, aa, freq_tot, freq_H, freq_E, freq_c):
	labels = aa
	
	x = np.arange(len(labels)) # the label locations
	width = 0.20  # the width of the bars

	fig, ax = plt.subplots()
	rects1 = ax.bar(x - width*3/2, freq_tot, width, label='overall')
	rects2 = ax.bar(x - width/2, freq_H, width, label='helix')
	rects3 = ax.bar(x + width/2, freq_E, width, label='strands')
	rects4 = ax.bar(x + width*3/2, freq_c, width, label='coil')

	# Add some text for labels, title and custom x-axis tick labels, etc.
	ax.set_ylabel('Frequency', fontsize=12)
	ax.set_title('Dssp Class Frequency By Residue', fontsize=12)
	ax.set_xticks(x)
	ax.set_xticklabels(labels)
	# custom font sizes
	ax.tick_params(axis='both', which='major', labelsize=12)
	plt.rcParams.update({'font.size': 12}) # influences legend
	plt.rc('axes', titlesize=12)
	plt.rc('axes', labelsize=12)
	ax.legend()

	#autolabel(rects1, ax)
	#autolabel(rects2, ax)
	#autolabel(rects3, ax)
	#autolabel(rects4, ax)

	fig.tight_layout()
	plt.savefig(outFile, dpi=None, facecolor='w', edgecolor='w',
		orientation='portrait', papertype=None, format='png',
		transparent=False, bbox_inches=None, pad_inches=0.1,
		metadata=None)
	plt.show()

def autolabel(rects, ax):
	"""Attach a text label above each bar in *rects*, displaying its height."""
	for rect in rects:
		height = rect.get_height()
		ax.annotate('{}'.format(height),
			xy=(rect.get_x() + rect.get_width() / 2, height),
			xytext=(0, 3),  # 3 points vertical offset
			textcoords="offset points",
			ha='center', va='bottom')
	
if __name__ == '__main__':
	sourceFile = sys.argv[1]
	outFile = sys.argv[2]
	aa, H, E, c = getCounts(sourceFile)
	aa, freq_tot, freq_H, freq_E, freq_c = getFrequencies(aa, H, E, c)
	multybarsHist(outFile, aa, freq_tot, freq_H, freq_E, freq_c)
