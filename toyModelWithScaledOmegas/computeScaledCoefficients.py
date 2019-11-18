#! /usr/bin/python3

import argparse
from math import sqrt, pi, fabs

helpString = """
"""

parser = argparse.ArgumentParser(description=helpString)
parser.add_argument('delta', metavar='d', type=float, help='interaction range minus diameter')
parser.add_argument('ecc', metavar='e', type=float, help='eccentricity')
parser.add_argument('eBB', metavar='eBB', type=float, help='')
parser.add_argument('eBS', metavar='eBS', type=float, help='')
parser.add_argument('eSS', metavar='eSS', type=float, help='')
##parser.add_argument('', metavar='', type=int, help='')
args = parser.parse_args()

ecc = args.ecc
delta = args.delta
eBB = args.eBB
eBS = args.eBS
eSS = args.eSS


HSradius = 0.5
HSdiameter = 1.0
bigRadius = HSradius + delta/2
patchRadius = bigRadius - ecc

def computeOmega(Ra, Rb, rab):
  # BKL paper, formula 18
  if ( rab > Ra+Rb ):
    return 0.
  elif ( rab <= fabs(Ra-Rb) ):
    return (1./HSradius**3)*min(Ra,Rb)**3;
  else:
    tempSum = (Ra**2-Rb**2)/(2.*rab)
    return (0.25/HSradius**3)*( (2.*Ra+tempSum+rab/2.)*(Ra-tempSum-rab/2.)**2 +
                (2.*Rb-tempSum+rab/2.)*(Rb+tempSum-rab/2.)**2 )



fBB = computeOmega(bigRadius, bigRadius, HSdiameter)
fBS = computeOmega(bigRadius, patchRadius, HSdiameter - ecc)
fSS = computeOmega(patchRadius, patchRadius, HSdiameter - 2*ecc)

print ( fBB, fBS, fSS)


cBB = eBB * fBB
cBS = eBS * fBS
cSS = eSS * fSS

print("scaled OUTPUT:\n")

print("{:.6f}".format(cBB).rjust(10) + "\t" + "{:.6f}".format(cBS).rjust(10) + "\t" +"{:.6f}".format(cBS).rjust(10) + "\n")
print("{:.6f}".format(cSS).rjust(10) + "\t" + "{:.6f}".format(cSS).rjust(10) + "\t" +"{:.6f}".format(cSS).rjust(10) + "\n")
#print(str(cBB).rjust(10) + "\t" + str(cBS).rjust(10) + "\t" + str(cBS).rjust(10) + "\n")
#print(str(cSS).rjust(10) + "\t" + str(cSS).rjust(10) + "\t" + str(cSS).rjust(10) + "\n")
print("1.00000\n")


cBB = eBB / fBB
cBS = eBS / fBS
cSS = eSS / fSS

print("multiplied OUTPUT:\n")
print(str(cBB).rjust(10) + "\t" + str(cBS).rjust(10) + "\t" + str(cBS).rjust(10) + "\n")
print(str(cSS).rjust(10) + "\t" + str(cSS).rjust(10) + "\t" + str(cSS).rjust(10) + "\n")
print("1.00000\n")
