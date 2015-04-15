#!/usr/bin/python

import numpy as np
import subprocess
import sys

def readFASTA(filename):
	headers = []
	strings = []
	with open(filename) as file:
		for line in file:
			if line.startswith('#'): continue
			if line.startswith('>'): 
				headers.append(line[1:len(line)-1])
				strings.append("")
			else: 
				strings[len(strings)-1] += line[:len(line)-1]
	return headers, strings

if len(sys.argv) == 1: print "No files to align"
for filename in sys.argv[1:]:
	headers, sequences = readFASTA(filename)
	#each sequence is a single cluster
	#calc first table
	A = np.zeros((len(sequences), len(sequences)))
	S = []
	for i in range(len(sequences)):
		for j in range(i + 1, len(sequences)):
			cmd = ['./pairwise', sequences[i], sequences[j]]
			process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
			result = []
			for line in process.stdout:
				#last symbol is \n
				result.append(line[:len(line)-1])
			if len(result) != 3:
				print "Error while trying to align sequences:"
				print sequences[i]
				print sequences[j]
				sys.exit()
			A[i, j] = int(result[2])
			S.append((result[0], result[1]))
	print A
	for elem in S:
		print elem
