#!/usr/bin/env python3

## TODO: https://pythonspot.com/matplotlib-pie-chart/
##	 https://matplotlib.org/3.1.1/gallery/pie_and_polar_charts/pie_features.html

import sys
import matplotlib.pyplot as plt

def readTable(fileName):
	x = []
	with open(fileName) as sourceFile:
		for line in sourceFile:
			x.append(line.rstrip())
	return(x)

def plotPieChart(data, labels, outFile, figSize):
	#fig = plt.figure(figsize=figSize)
	plt.pie(data, labels=labels, shadow=False, startangle=140)
	plt.axis('equal')
	plt.tight_layout()
	plt.show()
	plt.savefig(outFile, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)
	#return(fig)

def plotCircle(sizes, labels, outFile, figSize):
	colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99']
 
	fig1, ax1 = plt.subplots()
	ax1.pie(sizes, colors = colors, labels=labels, autopct='%1.1f%%', startangle=90)
	
	#draw circle
	centre_circle = plt.Circle((0,0),0.70,fc='white')
	fig = plt.gcf()
	fig.gca().add_artist(centre_circle)
	
	# Equal aspect ratio ensures that pie is drawn as a circle
	ax1.axis('equal')  
	plt.tight_layout()
	#plt.show()
	plt.savefig(outFile, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)

def donut(data, labels):
	import numpy as np
	fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))

	recipe = labels
	
	total = 0
	for x in data:
		total = total + float(x)
	percentages = []
	for i in range(len(data)):
		percentages.append(round(float(data[i]) / total * 100, 2))

	wedges, texts = ax.pie(data, wedgeprops=dict(width=0.5), startangle=0)
	
	for i in range(len(labels)):
		recipe[i] = labels[i] + '  (' + str(percentages[i]) + ' %)'

	#bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
	#kw = dict(arrowprops=dict(arrowstyle="-"),
	#	bbox=bbox_props, zorder=0, va="center")
	kw = dict(arrowprops=dict(arrowstyle="-"),zorder=0, va="center")

	for i, p in enumerate(wedges):
		ang = (p.theta2 - p.theta1)/2. + p.theta1
		y = np.sin(np.deg2rad(ang))
		x = np.cos(np.deg2rad(ang))
		horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
		connectionstyle = "angle,angleA=0,angleB={}".format(ang)
		kw["arrowprops"].update({"connectionstyle": connectionstyle})
		ax.annotate(recipe[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
			horizontalalignment=horizontalalignment, **kw)

	plt.tight_layout()
	plt.savefig(outFile, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)

if __name__ == '__main__':
	dataSource = sys.argv[1]
	labelSource = sys.argv[2]
	figSize = sys.argv[3]
	outFile = sys.argv[4]
	sizes = readTable(dataSource)
	labels = readTable(labelSource)
	#plotPieChart(data, labels, outFile, figSize)
	#plotCircle(sizes, labels, outFile, figSize)
	donut(sizes, labels)
