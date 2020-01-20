#!/usr/bin/env python3

import numpy as np
from sklearn.metrics import confusion_matrix
from sklearn.metrics import matthews_corrcoef
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def readDssp(dssp_dir, prot_id):
	dsspFile = open(dssp_dir + '/' + prot_id + '.dssp')
	dsspFile.readline()
	sequence = dsspFile.readline().rstrip()
	return(sequence)

#################
### confusion ###
#################

## MEMO:
## H = 0; E = 1; - = 2
def count_prediction(real_dssp_dir, imputed_dssp_dir, idList):
	real = []
	imputed = []
	with open(idList) as ids:
		for prot_id in ids:
			prot_id = prot_id.rstrip()
			real_dssp = readDssp(real_dssp_dir, prot_id)
			imputed_dssp = readDssp(imputed_dssp_dir, prot_id)
			for i in range(len(real_dssp)):
				if real_dssp[i] == 'H':
					real.append(0)
				if real_dssp[i] == 'E':
					real.append(1)
				if real_dssp[i] == '-':
					real.append(2)
				if imputed_dssp[i] == 'H':
					imputed.append(0)
				if imputed_dssp[i] == 'E':
					imputed.append(1)
				if imputed_dssp[i] == '-':
					imputed.append(2)
	real = np.array(real).astype(int)
	imputed = np.array(imputed).astype(int)
	return(real, imputed)

def build_conf_mat(y_true, y_pred):
	conf_mat = confusion_matrix(y_true, y_pred, labels=[0,1,2])
	return(conf_mat)

def print_confusion_matrix(confusion_matrix, class_names, figsize = (4,3), fontsize=14):
	df_cm = pd.DataFrame(
		confusion_matrix, index=class_names, columns=class_names, 
	)
	fig = plt.figure(figsize=figsize)
	try:
		heatmap = sns.heatmap(df_cm, annot=True, fmt="d")
	except ValueError:
		raise ValueError("Confusion matrix values must be integers.")
	heatmap.yaxis.set_ticklabels(heatmap.yaxis.get_ticklabels(), rotation=0, ha='right', fontsize=fontsize)
	heatmap.xaxis.set_ticklabels(heatmap.xaxis.get_ticklabels(), rotation=0, ha='right', fontsize=fontsize)
	heatmap.set_ylim(heatmap.get_ylim()[0] + 0.5, heatmap.get_ylim()[1] - 0.5)
	plt.ylabel('True label')
	plt.xlabel('Predicted label')
	plt.tight_layout()
	return(fig)

################
### identity ###
################

def compute_overall_identity(real_dssp_dir, imputed_dssp_dir, prot_id):
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

def compute_ss_identity(real_dssp_dir, imputed_dssp_dir, prot_id, ss):
	real = readDssp(real_dssp_dir, prot_id)
	imputed = readDssp(imputed_dssp_dir, prot_id)
	if len(real) != len(imputed):
		print('E: %s has length inconsistency'%(prot_id))
	else:
		identity = 0
		n_ss = 0
		for pos in range(len(real)):
			if real[pos] == ss:
				n_ss = n_ss + 1
				if real[pos] == imputed[pos]:
					identity = identity + 1
		if n_ss == 0:
			identity = -1
		else:
			identity = identity / n_ss
	return(identity)

def performance_identity(real_dssp_dir, imputed_dssp_dir, windowSize, idList):
	identities = {}
	with open(idList) as ids:
		for prot_id in ids:
			prot_id = prot_id.rstrip()
			identities[prot_id] = {}
			identity = compute_overall_identity(real_dssp_dir, imputed_dssp_dir, prot_id)
			identities[prot_id]['tot'] = identity
			for ss in ['H','E','-']:
				identities[prot_id][ss] = compute_ss_identity(real_dssp_dir, imputed_dssp_dir, prot_id, ss)
	return(identities)

def get_ss_identities(identities, ss):
	H = {}
	for prot_id in identities.keys():
		if identities[prot_id][ss] != -1:
			H[prot_id] = identities[prot_id][ss]
	return(H)

###########################################
### Overall Performance Across CV Folds ###

def overall_cv_folds_identity(real_dssp_dirs, imputed_dssp_dirs, idLists):
	real_dssp_dirs = real_dssp_dirs.split()
	imputed_dssp_dirs = imputed_dssp_dirs.split()
	idLists = idLists.split()
	overall_real = []
	overall_imputed = []
	for fold in range(len(idLists)):
		real_dssp_dir = real_dssp_dirs[fold]
		imputed_dssp_dir = imputed_dssp_dirs[fold]
		idList = idLists[fold]
		real, imputed = count_prediction(real_dssp_dir, imputed_dssp_dir, idList)
		for i in real:
			overall_real.append(i)
		for i in imputed:
			overall_imputed.append(i)
	return(overall_real, overall_imputed)

if __name__ == '__main__':
	import argparse as ap
	parser = ap.ArgumentParser(description='Compare the imputed dssp and the real dssp giving the gor performance')
	parser.add_argument('--imputed_dssp_dirs', dest='imputed_dssp_dirs', help='imputed dssp directories')
	parser.add_argument('--real_dssp_dirs', dest='real_dssp_dirs', help='real dssp directories')
	parser.add_argument('--id_lists', dest='id_lists', help='files containing lists of file names without extension')
	parser.add_argument('--out', dest='out_file', help='output file name')
	parser.add_argument('--window_size', dest='window_size', type=int, help='max number of residues that define the training windows size')
	args = parser.parse_args()
	imputed_dssp_dirs = args.imputed_dssp_dirs
	real_dssp_dirs = args.real_dssp_dirs
	idLists = args.id_lists
	outFile = args.out_file
	windowSize = args.window_size
	
	## test ##
	print(imputed_dssp_dirs)
	print()
	print(real_dssp_dirs)
	print()
	print(idLists)

	############################
	### identity performance ###
	
	## outdated
	#identities = performance_identity(real_dssp_dir, imputed_dssp_dir, windowSize, idList)
	#H = get_ss_identities(identities, 'H')
	#E = get_ss_identities(identities, 'E')
	#c = get_ss_identities(identities, '-')
	#print(np.mean(list(H.values())))
	#print(np.mean(list(E.values())))
	#print(np.mean(list(c.values())))

	########################
	### confusion matrix ###

	real, imputed = overall_cv_folds_identity(real_dssp_dirs, imputed_dssp_dirs, idLists)
	conf_mat = build_conf_mat(real, imputed)
	class_names = ['H','E','C']
	fig = print_confusion_matrix(conf_mat, class_names, fontsize=10)
	plt.savefig(outFile, type="png", dpi=96)

	###########
	### MCC ###

	print(matthews_corrcoef(real, imputed))
