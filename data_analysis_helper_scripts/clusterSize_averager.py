#! /usr/bin/python3

helpString = """Computes the line-by-line average of a specific file in the simulation directories.
This version was made extra for the cluster size file, which skips some of the sizes but only two columns. The results are read and collcted in a dictionary then printed out.
"""

import argparse
import os

parser = argparse.ArgumentParser(description=helpString)
parser.add_argument('file_name', metavar='f', type=str, help='file name')
args = parser.parse_args()

filename = args.file_name

dirs = [i for i in os.listdir() if i.startswith("siml_") ]

runs = [open(i + "/" + filename, "r") for i in dirs]

averaged = dict()

for run in runs:
  for fileline in run:
    a = fileline.split()
    if a[0][0] == "#":
      continue
    if( (a[0] in averaged) == False ):
      averaged[a[0]] = float(a[1])
    else:
      averaged[a[0]] += float(a[1])

outFile = open(filename + ".averaged", "w")
norm = 1./len(runs)
for i in averaged:
  outFile.write(str(i) + "\t" + str(norm*averaged[i]) + "\n")
