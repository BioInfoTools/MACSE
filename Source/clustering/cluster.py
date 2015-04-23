#!/usr/bin/python
# -*- coding: utf-8 -*- 

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

#alignment two blocks of sequences
#using clustal omega to do align
def align(data1, data2):
	#create subprocess
	child = subprocess.Popen(["clustalo","--in=-","--output-order=input-order"], \
		stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	#load sequences from first data block
	for i in range(len(data1)):
		if i % 2 == 0:
			child.stdin.write('>'+str(data1[i])+'\n')
		else:
			child.stdin.write(data1[i]+'\n')
	#load sequences from second data block
	for i in range(len(data2)):
		if i % 2 == 0:
			child.stdin.write('>'+str(data2[i])+'\n')
		else:
			child.stdin.write(data2[i]+'\n')
	#close stdin
	child.stdin.close()
	#save results
	index = -1
	result = data1 + data2
	for line in iter(child.stdout.readline,''):
		if line.startswith('>'): 
			index += 2
			result[index] = ""
		else: 
			if index > 0:
				result[index] += line.rstrip()
	return result

#m - distance between clusters
#s - align betwenn clusters
#n - count of clusters
def UPGMA(m, s, n):
	if n == 2:
		return s[1][0] # <- answer
	#search for closest clusters
	max_value, index1, index2 = 0, 0, 0
	for i in range(n):
		for j in range(i):
			if m[i][j] > max_value:
				index1, index2, max_value = i, j, m[i][j]
	tek_align = s[index1][index2]
	#updating distance & aligns
	for i in range(n):
		for j in range(i):
			if not (i == index1 or j == index1) and (i == index2 or j == index2):
				m[i][j] = (m[i][j] + max_value) / 2.0
				s[i][j] = align(s[i][j], tek_align)
	#delete columns & raws
	for i in range(index1+1, n):
		del m[i][index1]
		del s[i][index1]
	del m[index1]
	del s[index1]
	#next iteration
	return UPGMA(m, s, n-1)

#MAIN
if len(sys.argv) == 1: print "No files to align"
for filename in sys.argv[1:]:
	#parsing input file
	headers, sequences = readFASTA(filename)
	#creating distance matrix
	#each sequence is a single cluster
	matrix = []
	aligns = []
	for i in range(0, len(sequences)):
		matrix.append([])
		aligns.append([])
		for j in range(0, i):
			cmd = ['./pairwise', sequences[i], sequences[j]]
			process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
			result = []
			for line in process.stdout:
				result.append(line.rstrip())
			#check for an error
			if len(result) != 3:
				print "Error while trying to align sequences:"
				print headers[i]
				print headers[j]
				print ''.join(map(lambda x: x + ' ', cmd))
				sys.exit()
			#save data
			aligns[i].append([headers[i], result[0], headers[j], result[1]])
			matrix[i].append(int(result[2]))
		matrix[i].append(0)
	#clustering
	UPGMA(matrix, aligns, len(sequences))
	