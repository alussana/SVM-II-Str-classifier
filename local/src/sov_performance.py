#!/opt/conda/bin/python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import seaborn as sns
from sklearn.metrics import accuracy_score

def readDssp(dssp_dir, prot_id):
	dsspFile = open(dssp_dir + '/' + prot_id + '.dssp')
	dsspFile.readline()
	sequence = dsspFile.readline().rstrip()
	return(sequence)

def compute_identity(real_dssp_dir, imputed_dssp_dir, prot_id):
	real = readDssp(real_dssp_dir, prot_id)
	imputed = readDssp(imputed_dssp_dir, prot_id)
	if len(real) != len(imputed):
		print('E: %s has length inconsistency'%(prot_id))
	else:
		identity = 0
		for pos in range(len(real)):
			if real[pos] == imputed[pos]:
				identity = identity + 1
		identity = identity / len(real)
	return(identity)

def performance_identity(real_dssp_dir, imputed_dssp_dir, windowSize, idList):
	identities = {}
	with open(idList) as ids:
		for prot_id in ids:
			prot_id = prot_id.rstrip()
			identity = compute_identity(real_dssp_dir, imputed_dssp_dir, prot_id)
			identities[prot_id] = identity
	return(identities)

def print_identities(identities, outFile):
	out = open(outFile, 'w')
	Q3 = []
	for prot_id in identities:
		print('%s\t%f'%(prot_id, identities[prot_id]), file=out)
		Q3.append(identities[prot_id])
	mean_Q3 = np.mean(Q3)
	std_err_Q3 = np.std(Q3) / np.sqrt(len(Q3))
	print(f'mean_Q3: {mean_Q3}\tstd_err_Q3: {std_err_Q3}')

def compute_sov_segments(real, imputed, ss):
	real_segments = []
	imputed_segments = []
	pos = 0; real_start = -1; imputed_start = -1; real_end = -1; imputed_end = -1
	while pos < len(real):
		if real[pos] == ss:
			real_end = pos
			if real_start == -1:
				real_start = pos
		else:
			if real_end != -1:
				real_segments.append((real_start, real_end))
				real_start = -1
				real_end = -1
		if imputed[pos] == ss:
			imputed_end = pos
			if imputed_start == -1:
				imputed_start = pos 
		else:
			if imputed_end != -1:
				imputed_segments.append((imputed_start, imputed_end))
				imputed_start = -1
				imputed_end = -1
		pos = pos + 1
	if real_end != -1:
		real_segments.append((real_start, real_end))
	if imputed_end != -1:
		imputed_segments.append((imputed_start, imputed_end))
	return(real_segments, imputed_segments)

def compute_sov(real_dssp_dir, imputed_dssp_dir, prot_id):
	real = readDssp(real_dssp_dir, prot_id)
	imputed = readDssp(imputed_dssp_dir, prot_id)
	if len(real) != len(imputed):
		print('E: %s has length inconsistency'%(prot_id))
	else:
		sov = {'H': 0,'E': 0,'-': 0}
		N = {'H': 0,'E': 0,'-': 0}
		for ss in ['H','E','-']:
			## divide sequences into segments
			real_segments, imputed_segments = compute_sov_segments(real, imputed, ss)
			## if no segments in conformation ss are found in reference dssp, then continue
			if len(real_segments) == 0: 
				continue
			## compute sov[ss]
			for real_seg in real_segments:
				## check overlap with any imputed segment
				## compute minov and maxov
				minov = 0
				maxov = 0
				N_updated = 0
				for imputed_seg in imputed_segments:
					## determine intersection of residues
					real_pos = set(range(real_seg[0], real_seg[1] + 1))
					imputed_pos = set(range(imputed_seg[0], imputed_seg[1] + 1))
					intersection = real_pos & imputed_pos
					## if intersection is not null then compute minov, maxov, len_s1, delta
					if len(intersection) > 0:
						minov = len(intersection)
						maxov = len(real_pos) + len(imputed_pos) - len(intersection)
						len_s1 = len(set(range(real_seg[0], real_seg[1] + 1)))
						delta = min(maxov - minov, minov, int(len(real_pos) / 2), int(len(imputed_pos) / 2))
						## update sov[ss]
						sov[ss] = sov[ss] + len_s1 * (minov + delta) / maxov
						## update N[ss]
						N[ss] = N[ss] + len_s1
						N_updated = 1
						## test ##
						#print(f'{prot_id}\tss: {ss}\tminov: {minov}\tmaxov: {maxov}\tdelta: {delta}\tN: {N[ss]}')
						#########
				## update N[ss] if real segment has no overlaps
				if N_updated == 0:
					N[ss] = N[ss] + len(set(range(real_seg[0], real_seg[1] + 1)))
			## test ##
			#print(f'{prot_id}\tss: {ss}\tsov: {sov[ss]}')
			##########
		## aggregate sov score for all the secondary structures
		total_N = N['H'] + N['E'] + N['-']
		total_sov = sov['H'] + sov['E'] + sov['-']
		sov99 = 100 * total_sov / total_N
		for ss in sov.keys():
			if N[ss] == 0:
				sov[ss] = 'NA'
			else:
				sov[ss] = 100 * sov[ss] / N[ss]
		return(sov99, sov)

def performance_sov(real_dssp_dir, imputed_dssp_dir, windowSize, idList):
	sov = {}
	with open(idList) as ids:
		for prot_id in ids:
			prot_id = prot_id.rstrip()
			sov_99, partial_sov = compute_sov(real_dssp_dir, imputed_dssp_dir, prot_id)
			sov[prot_id] = {}
			for ss in partial_sov.keys():
				sov[prot_id][ss] = partial_sov[ss]
			sov[prot_id]['99_3'] = sov_99
	return(sov)

def print_sov(sov, outFile):
	out = open(outFile, 'w')
	scores = {}
	mean_sov = {}
	std_err_sov = {}
	for score in sov[list(sov.keys())[0]]:
		scores[score] = []
		mean_sov[score] = 0
		std_err_sov[score] = 0
	for prot_id in sov.keys():
		toPrint = prot_id
		for score in sov[prot_id]:
			toPrint = toPrint + f'\t{score}: {sov[prot_id][score]}'
			if sov[prot_id][score] != 'NA':
				scores[score].append(sov[prot_id][score])
		print(toPrint, file = out)
	out.close()
	for score in scores:
		mean_sov[score] = np.mean(scores[score])
		std_err_sov[score] = np.std(scores[score]) / np.sqrt(len(scores[score]))
	print('SOV')
	for score in mean_sov:
		m = mean_sov[score]
		e = std_err_sov[score]
		if score == '-':
			score = 'C'
		print(f'\t{score}: {m}\t{e}')

def feed_confusion_matrix(real_dssp_dir, imputed_dssp_dir, idList):
	y_real = []
	y_pred = []
	with open(idList) as ids:
		for prot_id in ids:
			prot_id = prot_id.rstrip()
			dssp_real = readDssp(real_dssp_dir, prot_id)
			dssp_pred = readDssp(imputed_dssp_dir, prot_id)
			if len(dssp_real) != len(dssp_pred):
				print('%s\tError: dssp length inconsistency'%(prot_id))
			for i in range(len(dssp_real)):
				if dssp_real[i] == 'H':
					y_real.append(0)
				if dssp_real[i] == 'E':
					y_real.append(1)
				if dssp_real[i] == '-':
					y_real.append(2)
			for i in range(len(dssp_pred)):
				if dssp_pred[i] == 'H':
					y_pred.append(0)
				if dssp_pred[i] == 'E':
					y_pred.append(1)
				if dssp_pred[i] == '-':
					y_pred.append(2)
	return(y_real, y_pred)

#def print_confusion_matrix(confusion_matrix, class_names, figsize = (4,3), fontsize=14):
#	"""Prints a confusion matrix, as returned by sklearn.metrics.confusion_matrix, as a heatmap.
#    	
#	Arguments
#	---------
#	confusion_matrix: numpy.ndarray
#	    The numpy.ndarray object returned from a call to sklearn.metrics.confusion_matrix. 
#	    Similarly constructed ndarrays can also be used.
#	class_names: list
#	    An ordered list of class names, in the order they index the given confusion matrix.
#	figsize: tuple
#	    A 2-long tuple, the first value determining the horizontal size of the ouputted figure,
#	    the second determining the vertical size. Defaults to (10,7).
#	fontsize: int
#	    Font size for axes labels. Defaults to 14.
#	    
#	Returns
#	-------
#	matplotlib.figure.Figure
#	    The resulting confusion matrix figure
#	"""
#	df_cm = pd.DataFrame(
#		confusion_matrix, index=class_names, columns=class_names, 
#	)
#	fig = plt.figure(figsize=figsize)
#	try:
#		heatmap = sns.heatmap(df_cm, annot=True, fmt="d")
#	except ValueError:
#		raise ValueError("Confusion matrix values must be integers.")
#	heatmap.yaxis.set_ticklabels(heatmap.yaxis.get_ticklabels(), rotation=0, ha='right', fontsize=fontsize)
#	heatmap.xaxis.set_ticklabels(heatmap.xaxis.get_ticklabels(), rotation=0, ha='right', fontsize=fontsize)
#	plt.ylabel('True label')
#	plt.xlabel('Predicted label')
#	plt.tight_layout()
#	return(fig)

def print_blue_confusion_matrix(y_true, y_pred, classes, normalize=True, title=None, cmap=plt.cm.Blues):
	from sklearn.utils.multiclass import unique_labels
	cm = confusion_matrix(y_true, y_pred)
	#classes = classes[unique_labels(y_true, y_pred)]
	#classes = [0, 1, 2]
	classes = ['H','E','C']
	if normalize:
		cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
	print(cm)
	fig, ax = plt.subplots()
	im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
	ax.figure.colorbar(im, ax=ax).ax.tick_params(labelsize=12)
	## custom tick_params (test)
	ax.tick_params(axis='both', which='major', labelsize=16)
	plt.rcParams.update({'font.size': 16})
	#plt.rc('font', size=18)	# controls default text sizes
	#plt.rc('axes', titlesize=18)	# fontsize of the axes title
	#plt.rc('axes', labelsize=18)	# fontsize of the x and y labels
	#plt.rc('xtick', labelsize=18)	# fontsize of the tick labels
	#plt.rc('ytick', labelsize=18)	# fontsize of the tick labels
	#plt.rc('legend', fontsize=18)	# legend fontsize
	#plt.rc('figure', titlesize=18)	# fontsize of the figure title
	# We want to show all ticks...
	ax.set(xticks=np.arange(cm.shape[1]),
		yticks=np.arange(cm.shape[0]),
		# ... and label them with the respective list entries
		xticklabels=classes, yticklabels=classes,
		title=title,
		ylabel='True label',
		xlabel='Predicted label')
	## custom fix for the matplotlib bug (https://github.com/matplotlib/matplotlib/issues/14751)
	ax.set_ylim(len(classes)-0.5, -0.5)
	## custom fix: override the labels using the wanted font size
	plt.xlabel('Predicted Class', fontsize=16)
	plt.ylabel('Observed Class', fontsize=16)
	# Rotate the tick labels and set their alignment.
	plt.setp(ax.get_xticklabels(), rotation=0, ha="center",
		rotation_mode="anchor")
	# Loop over data dimensions and create text annotations.
	fmt = '.2f' if normalize else 'd'
	thresh = cm.max() / 2.
	for i in range(cm.shape[0]):
		for j in range(cm.shape[1]):
			ax.text(j, i, format(cm[i, j], fmt),
				ha="center", va="center",
				color="white" if cm[i, j] > thresh else "black")
	fig.tight_layout()
	return ax

def analyze_conf_mat(confMat, class_dict):
	analysis = {}
	for ss in class_dict:
		analysis[ss] = {}
		index = class_dict[ss]
		non_index = list(class_dict.values())
		non_index.pop(index)
		analysis[ss]['TP'] = confMat[index][index]
		TN = 0
		for x in non_index:
			for y in non_index:
				TN = TN + confMat[x][y]
		analysis[ss]['TN'] = TN
		FP = 0
		for x in non_index:
			FP = FP + confMat[x][index]
		analysis[ss]['FP'] = FP
		FN = 0
		for x in non_index:
			FN = FN + confMat[index][x]	
		analysis[ss]['FN'] = FN
	return(analysis)

def compute_indexes(analysis):
	indexes = {}
	for ss in analysis:
		indexes[ss] = {}
		TP = analysis[ss]['TP']
		TN = analysis[ss]['TN']
		FP = analysis[ss]['FP']
		FN = analysis[ss]['FN']
		indexes[ss]['TPR'] = TP / (TP + FN)
		indexes[ss]['TNR'] = TN / (TN + FP)
		indexes[ss]['PPV'] = TP / (TP + FP)
		indexes[ss]['NPV'] = TN / (TN + FN)
		indexes[ss]['FPR'] = FP / (FP + TN)
		indexes[ss]['FNR'] = FN / (TP + FN)
		indexes[ss]['FDR'] = FP / (TP + FP)
		indexes[ss]['MCC'] = (TP * TN - FP * FN) / (np.sqrt((TP + FP)* (TP + FN) * (TN + FP) * (TN + FN)))
	return(indexes)

def print_indexes(indexes):
	for ss in indexes:
		MCC = indexes[ss]['MCC']
		PPV = indexes[ss]['PPV']
		TPR = indexes[ss]['TPR']
		print(ss)
		print(f'\tMCC: {MCC}')
		print(f'\tPPV: {PPV}')
		print(f'\tTPR: {TPR}')

if __name__ == '__main__':
	import argparse as ap
	parser = ap.ArgumentParser(description='Compare the imputed dssp and the real dssp giving the gor performance')
	parser.add_argument('--imputed_dssp_dir', dest='imputed_dssp_dir', help='imputed dssp directory')
	parser.add_argument('--real_dssp_dir', dest='real_dssp_dir', help='real dssp directory')
	parser.add_argument('--id_list', dest='id_list', help='list of file names without extension')
	parser.add_argument('--out', dest='out_file', help='output file name')
	parser.add_argument('--window_size', dest='window_size', type=int, help='max number of residues that define the training windows size')
	args = parser.parse_args()
	imputed_dssp_dir = args.imputed_dssp_dir
	real_dssp_dir = args.real_dssp_dir
	idList = args.id_list
	outFile = args.out_file
	windowSize = args.window_size

	############################
	### identity performance ###
	identities = performance_identity(real_dssp_dir, imputed_dssp_dir, windowSize, idList)
	print_identities(identities, outFile)
	
	########################
	### Confusion Matrix ###
	class_names = ['H','E','C']
	y_true, y_pred = feed_confusion_matrix(real_dssp_dir, imputed_dssp_dir, idList)
	fig = print_blue_confusion_matrix(y_true, y_pred, class_names)
	confMat = confusion_matrix(y_true, y_pred)
	figname = imputed_dssp_dir + '_confusion_matrix.png'
	plt.savefig(figname, type="png", dpi=96)
	print(f'acc_score: {accuracy_score(y_true, y_pred)}')
	
	###############
	### Indexes ###
	pointers = [0, 1, 2]
	class_dict = dict(zip(class_names, pointers))
	analysis = analyze_conf_mat(confMat, class_dict)
	indexes = compute_indexes(analysis)
	print_indexes(indexes)

	#######################
	### sov performance ###
	sov = performance_sov(real_dssp_dir, imputed_dssp_dir, windowSize, idList)
	print_sov(sov, outFile)
