#!/usr/bin/env python3

import sys
import gzip
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def readTable(fileName):
	#domainID = []
	sequence = []
	with gzip.open(fileName, 'rt') as sourceTable:
		for sourceLine in sourceTable:
			sourceInfo = sourceLine.split("\t")
			#domainID.append(sourceInfo[0])
			sequence.append(sourceInfo[1])
	return(sequence)

def computePropensity(windowSize, aa, fasta, dssp):
	H = []
	E = []
	for i in range(windowSize):
		H.append(0)
		E.append(0)
	for x in range(len(fasta)):
		i = 0
		n = len(fasta[x]) - 1
		while n - i >= windowSize:
			if dssp[x][i + windowSize // 2] == 'H':
				for j in range(windowSize):
					if fasta[x][i+j] == aa:
						H[j] = H[j] + 1
			elif dssp[x][i + windowSize // 2] == 'E':
				for j in range(windowSize):
					if fasta[x][i+j] == aa:
						E[j] = E[j] + 1
			i = i + 1
	return((H,E))

def generateTables(multitable):
	helixTable = []
	strandTable = []
	for aa in multitable:
		helixTable.append(aa[0])
		strandTable.append(aa[1])
	return(np.array(helixTable).astype(float), np.array(strandTable).astype(float))

def calculateHeatmap(fastaSource, dsspSource, outFile, windowSize):
	aa = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']
	fasta = readTable(fastaSource)
	dssp = readTable(dsspSource)
	propensity17 = []
	for aminoAcid in aa:
		propensity17.append(computePropensity(windowSize, aminoAcid, fasta, dssp))
	helixTable, strandTable = generateTables(propensity17)
	tot_counts = 0
	for j in range(len(aa)):
		tot_counts = tot_counts + helixTable[j][0]
	for i in range(windowSize):
		for j in range(len(aa)):
			helixTable[j][i] = helixTable[j][i] / tot_counts
			strandTable[j][i] = strandTable[j][i] / tot_counts
	return(helixTable, strandTable)

'''
def doubleHeatmap(aa, propensityTable, title, outFile):
	vegetables = aa
	farmers = [i for i in range(17)]
	harvest = propensityTable
	fig, ax = plt.subplots()
	im = ax.imshow(harvest)

	# show all ticks
	ax.set_xticks(np.arange(len(farmers)))
	ax.set_yticks(np.arange(len(vegetables)))

	# label them with the respective list entries
	ax.set_xticklabels(farmers)
	ax.set_yticklabels(vegetables)

	# rotate the tick labels and set their alignment
	plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

	# Create colorbar
	cbar = ax.figure.colorbar(im, ax=ax)
	cbar.ax.set_ylabel(frequency, rotation=-90, va="bottom")

	# Loop over data dimensions and create text annotations.
	#for i in range(len(vegetables)):
    	#	for j in range(len(farmers)):
        #		text = ax.text(j, i, harvest[i][j], ha="center", va="center", color="w")

	ax.set_title(title)
	fig.tight_layout()
	plt.show()
	plt.savefig(outFile, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)
'''

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    #cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    #cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    #ax.tick_params(top=True, bottom=False,
    #               labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), ha="center",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    #return im, cbar
    return im


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

if __name__ == '__main__':
	# fastaSource and dsspSource have to be consistently sorted
	#aa = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']
	aa = ['G', 'A', 'V', 'P', 'L', 'I', 'M', 'F', 'W', 'Y', 'S', 'T', 'C', 'N', 'Q', 'H', 'D', 'E', 'K', 'R']
	fastaSource = sys.argv[1]
	dsspSource = sys.argv[2]
	outFile = sys.argv[3]
	windowSize = 17
	helixTable, strandTable = calculateHeatmap(fastaSource, dsspSource, outFile, windowSize)

	## PLOT ##
	if outFile.endswith('helix.png'):
		title = 'Residue Composition: Helix'
		propensityTable = helixTable
	elif outFile.endswith('strand.png'):
		title = 'Residue Composition: Strand'
		propensityTable = strandTable
	else:
		title = 'Residue Composition'
		propensityTable = helixTable
	#doubleHeatmap(aa, propensityTable, title, outFile)
	fig, ax = plt.subplots()
	farmers = [(i-windowSize // 2) for i in range(len(aa))]
	#im, cbar = heatmap(propensityTable, aa, farmers, ax=ax,
        #	           cmap="viridis", cbarlabel="Frequency")
	im = heatmap(propensityTable, aa, farmers, ax=ax, cmap="viridis", cbarlabel="Frequency")
	#texts = annotate_heatmap(im, valfmt="{x:.1f} t")
	## trying to set the font size
	ax.figure.colorbar(im, ax=ax).ax.tick_params(labelsize=12)
	ax.tick_params(axis='both', which='major', labelsize=12)
	fig.tight_layout()
	#plt.show()
	plt.savefig(outFile, dpi=None, facecolor='w', edgecolor='w',
		orientation='portrait', papertype=None, format='png',
		transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)
