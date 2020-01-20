#! /usr/bin/python3

helpString = """Computes the line-by-line average of a specific file in the simulation directory.
"""

import argparse
import os

parser = argparse.ArgumentParser(description=helpString)
parser.add_argument('file_name', metavar='f', type=str, help='file name')
parser.add_argument('columns', metavar='cN', type=int, help='number of columns in file')
args = parser.parse_args()

filename = args.file_name
numColumns = args.columns


dirs = [i for i in os.listdir() if i.startswith("siml_") ]

files = [open(i + "/" + filename, "r") for i in dirs]

averaged = []

for rows in zip(*files):
  averagedLine = [0.0 for i in range(numColumns)]
  for fileline in rows:
    a = fileline.split()
    if a[0][0] == "#":
      continue
    for i in range(len(a)):
      averagedLine[i] += float(a[i])
  for i in range(len(averagedLine)):
    averagedLine[i] /= len(dirs)
  averaged.append(averagedLine)

outFile = open(filename + ".averaged", "w")
for i in averaged:
  for j in i:
    outFile.write(str(j) + "\t")
  outFile.write("\n")
