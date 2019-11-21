#! /usr/bin/python3

import argparse
from math import sqrt, pi, fabs

helpString = """
"""

parser = argparse.ArgumentParser(description=helpString)
parser.add_argument('mode', metavar='m', type=str, help='usage mode')
parser.add_argument('delta', metavar='d', type=float, help='interaction range minus diameter')
parser.add_argument('ecc', metavar='e', type=float, help='eccentricity')
parser.add_argument('eBB', metavar='eBB', type=float, help='')
parser.add_argument('eBS', metavar='eBS', type=float, help='')
parser.add_argument('eSS', metavar='eSS', type=float, help='')
parser.add_argument('emin', metavar='emin', type=float, help='')
##parser.add_argument('', metavar='', type=int, help='')
args = parser.parse_args()

usageMode = args.mode
ecc = args.ecc
delta = args.delta
eBB = args.eBB
eBS = args.eBS
eSS = args.eSS
emin = args.emin

HSradius = 0.5
HSdiameter = 1.0
bigRadius = HSradius + delta/2
patchRadius = bigRadius - ecc

def computeOmega(Ra, Rb, rab):
  # BKL paper, formula 18
  if ( rab > Ra+Rb ):
    return 0.
  elif ( rab <= fabs(Ra-Rb) ):
    return (min(Ra,Rb)/HSradius)**3;
  else:
    tempSum = (Ra**2-Rb**2)/(2.*rab)
    return (0.25/HSradius**3) * (
               (2.*Ra+tempSum+rab/2.)*(Ra-tempSum-rab/2.)**2 +
               (2.*Rb-tempSum+rab/2.)*(Rb+tempSum-rab/2.)**2 )


fBB = computeOmega(bigRadius, bigRadius, HSdiameter)
fBS = computeOmega(bigRadius, patchRadius, HSdiameter - ecc)
fSS = computeOmega(patchRadius, patchRadius, HSdiameter - 2*ecc)

print("volumes at contact:")
print( fBB, fBS, fSS)

if usageMode == "include":
  cBB = eBB / (fBB * emin)
  cBS = eBS / (fBS * emin)
  cSS = eSS / (fSS * emin)
  print("scaled OUTPUT:")
elif usageMode == "extract":
  cBB = eBB * fBB / emin
  cBS = eBS * fBS / emin
  cSS = eSS * fSS / emin
  print("multiplied OUTPUT:")
else:
  print("mode does not exist")
  exit()


print(str(cBB)[0:7].ljust(10) + str(cBS)[0:8].ljust(10) + str(cBS)[0:8].ljust(10))
print(str(cSS)[0:7].ljust(10) + str(cSS)[0:7].ljust(10) + str(cSS)[0:7].ljust(10))

#print('{0: 1.5e}    {1: 1.5e}    {1: 1.5e}'.format(cBB, cBS))
#print('{0: 1.5e}    {0: 1.5e}    {0: 1.5e}'.format(cSS))
