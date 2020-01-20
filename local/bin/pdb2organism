#!/usr/bin/env python3

import sys
import requests
import xml.dom.minidom as minidom

def readTable(fileName):
	x = []
	with open(fileName) as sourceFile:
		for line in sourceFile:
			x.append(line.rstrip())
	return(x)

def getOrganism(IDs):
	for proteinId in IDs:
		urlForReport = 'http://www.rcsb.org/pdb/rest/customReport.xml?pdbids='
		result = requests.get(urlForReport + proteinId + "&customReportColumns=taxonomy")
		xmldoc = minidom.parseString(result.text)
		element = xmldoc.getElementsByTagName("dimEntity.taxonomy")
		print(element[0].childNodes[0].data)

if __name__ == '__main__':
	idsSource = sys.argv[1]
	IDs = readTable(idsSource)
	getOrganism(IDs)
