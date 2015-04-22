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

def UPGMA(A, N, S, i = 0, j = 0):
	def Align(s1, s2, subst="BLOSUM62"):
		#reading score matrix
		score_matrix = np.zeros((256, 256))
		with open(substr) as file:
			for line in file:
				if line.startswith('#'): continue
				#reading alphabet
				alphabet = line.replace(' ', '')
				break
			#reading scores
			for line in file:
				
		#initialize
		F = np.zeros((len(s1[0])+1, len(s2[0])+1))
		W = np.zeros((len(s1[0])+1, len(s2[0])+1))
		#forward pass
		for i in range(1, len(s1[0])+1):
			for j in range(1, len(s2[0])+1):
				way1 = F[i-1, j-1] + score_matrix[seq1[i-1]*128 + seq2[j-1]];
			//порвать seq2
			int way2 = F[(i-1)*m + j] + gap_extension;
			if (W[(i-1)*m + j] == 1) way2 += gap_open;
			//порвать seq1
			int way3 = F[i*m + j-1] + gap_extension;
			if (W[i*m + j-1] == 1) way3 += gap_open;
			//выбираем максимум
			if (way1 >= way2 && way1 >= way3) {
				F[i*m + j] = way1;
				W[i*m + j] = 1;
			} else if (way2 >= way1 && way2 >= way3) {
				F[i*m + j] = way2;
				W[i*m + j] = 2;
			} else {
				F[i*m + j] = way3;
				W[i*m + j] = 3;
			}
	if N == 1:
		#finish
		return S[0]
	if i == j:
		#search center clustering
		max_value = A[0, 1]
		for index1 in range(0, N):
			for index2 in range(i + 1, N):
				if max_value < A[index1, index2]:
					i = index1
					j = index2
					max_value = A[i, j]
	#creating a new cluster (i, j)
	cls = Align(s[i], s[j])

#MAIN
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
				print ''.join(map(lambda x: x + ' ', cmd))
				sys.exit()
			A[i, j] = int(result[2])
			S.append((result[0], result[1]))
			#Debug print
			#print ''.join(map(lambda x: x + ' ', cmd))
	
	#Debug out
	#print A
	#for elem in S:
	#	print elem
