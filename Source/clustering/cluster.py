#!/usr/bin/python

import numpy as np
import subprocess
import sys

def read_FASTA_strings(filename):
	with open(filename) as file:
		return file.read().split('>')[1:]

def read_FASTA_sequences(filename):
	return map(lambda seq: ''.join(seq.split('\n')[1:]), read_FASTA_strings(filename))

for filename in sys.argv[1:]:
	sequences = read_FASTA_sequences(filename)
	#each sequence is a single cluster
	#calc first table
	A = np.zeros((len(sequences), len(sequences)))
	for i in range(len(sequences)):
		for j in range(i + 1, len(sequences)):
			cmd = ['./pairwise', sequences[i], sequences[j]]
			process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
			result = []
			for line in process.stdout:
				#last symbol is \n
				result.append(line[:len(line)-1])
			print result[0]
			print result[1]
			A[i, j] = int(result[2])
	print A

