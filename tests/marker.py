#!/usr/bin/python
# -*- coding: utf-8 -*- 

import string
import sys

file1 = sys.argv[1]
file2 = sys.argv[2]

def ReadData(filename):
  with open(filename) as f:
    index = 0
    a = {}
    content = map(lambda x: x.rstrip(), f.readlines())
    while index < len(content):
      if content[index] == "":
        index += 1
        continue
      a[content[index]] = content[index+1]
      index += 2
  return a

data1 = ReadData(file1)
data2 = ReadData(file2)

total = 0.
for name in data1.keys():
  index1 = 0
  index2 = 0
  score = 0.
  seq = data1[name]
  while index1 < len(seq) and index2 < len(data2[name]):
    #print index1, index2, seq[index1], data2[name][index2]
    if (seq[index1] == data2[name][index2]):
      score += 1
      index1 += 1
      index2 += 1
    else:
      if (seq[index1] == '-'):
        index1 += 1
      else:
        if (data2[name][index2] == '-'):
          index2 += 1
        else:
          index1 += 1
          index2 += 1
  #print (score / len(seq) + score / len(data2[name])) / 2.
  total += (score / len(seq) + score / len(data2[name])) / 2.
  
print total / len(data1)
