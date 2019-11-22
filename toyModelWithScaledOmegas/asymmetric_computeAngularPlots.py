#! /usr/bin/python3

import argparse
from math import sqrt, pi, fabs, cos
from math import sin as sen

helpString = """Usage Modes: yes, no
Example: no 0.2 0.22 .245728   -3.11694  21.2298   .142495 (45n)
no 0.2 0.22  3.18802   -24.3562  58.9717   1.02726 (45c)
yes 0.2 0.22 -1 0 0 1
yes 0.2 0.22 0 0 -1 1
yes 0.2 0.22 0 -1 2 1
"""

parser = argparse.ArgumentParser(description=helpString)
parser.add_argument('mode', metavar='m', type=str, help='do the e_ij include the omega-scaling?')
parser.add_argument('delta', metavar='d', type=float, help='interaction range minus diameter')
parser.add_argument('ecc1', metavar='e1', type=float, help='eccentricity first patch')
parser.add_argument('ecc2', metavar='e2', type=float, help='eccentricity second patch')
parser.add_argument('eBB', metavar='eBB', type=float, help='')
parser.add_argument('eBs1', metavar='eBs1', type=float, help='')
parser.add_argument('eBs2', metavar='eBs2', type=float, help='')
parser.add_argument('es1s1', metavar='es1s1', type=float, help='')
parser.add_argument('es1s2', metavar='es1s2', type=float, help='')
parser.add_argument('es2s2', metavar='es2s2', type=float, help='')
parser.add_argument('emin', metavar='emin', type=float, help='')
##parser.add_argument('', metavar='', type=int, help='')
args = parser.parse_args()

usageMode = args.mode
ecc = args.ecc
delta = args.delta
eBB = args.eBB
eBs1 = args.eBs1
eBs2 = args.eBs2
es1s1 = args.es1s1
es1s2 = args.es1s2
es2s2 = args.es2s2
emin = args.emin

HSradius = 0.5
HSdiameter = 1.0
bigRadius = HSradius + delta/2
patch1radius = bigRadius - ecc1
patch2radius = bigRadius - ecc2

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

def computePotentialVerticalParticle(theta):
  rP1P2 = sqrt( (HSdiameter + ecc1*sen(theta))**2 + (ecc1*(1 - cos(theta)))**2 )
  rP1B2 = sqrt( ecc1*2 + HSdiameter**2 )
  rP1Q2 = sqrt( (HSdiameter - ecc2*sen(theta))**2 + (ecc1 + ecc2*cos(theta))**2 )
  rB1P2 = sqrt( (HSdiameter + ecc1*sen(theta))**2 + (ecc1*cos(theta))**2 )
  rB1B2 = HSdiameter
  rB1Q2 = sqrt( (HSdiameter - ecc2*sen(theta))**2 + (ecc2*cos(theta))**2 )
  rQ1P2 = sqrt( (HSdiameter + ecc1*sen(theta))**2 + (ecc2 + ecc1cos(theta))**2 )
  rQ1B2 = sqrt( ecc2*2 + HSdiameter**2 )
  rQ1Q2 = sqrt( (HSdiameter - ecc2*sen(theta))**2 + (ecc2*(1 - cos(theta)))**2 )

  wBB  = computeOmega(bigRadius, bigRadius,    rB1B2)
  wBs1 = computeOmega(bigRadius, patch1Radius, rB1Q2)
  wBs2 = computeOmega(bigRadius, patch2Radius, rB1Q2)
  
  V = (
        eBB*    computeOmega(bigRadius,    bigRadius,    rB1B2) +
        eBs1*(  computeOmega(patch1Radius, bigRadius,    rP1B2) +
                computeOmega(bigRadius,    patch1Radius, rB1P2) ) +
        eBs2*(  computeOmega(bigRadius,    patch2Radius, rB1Q2) +
                computeOmega(patch2Radius, bigRadius,    rQ1B2) ) +
        es1s1*  computeOmega(patch1Radius, patch1Radius, rP1P2) +
        es1s2*( computeOmega(patch1Radius, patch2Radius, rP1Q2) +
                computeOmega(patch2Radius, patch1Radius, rQ1P2) ) +
        es2s2*  computeOmega(patch2Radius, patch2Radius, rQ1Q2)
      )
  return V, wBB, wBs1

def computePotentialHorizontalParticle(theta):
  rP1P2 = sqrt( (HSdiameter + ecc1*(1 + sen(theta)))**2 + (ecc1*cos(theta))**2 )
  rP1B2 = ecc1 + HSdiameter
  rP1Q2 = sqrt( (HSdiameter - ecc1 - ecc2*sen(theta))**2 + (ecc2*cos(theta))**2 )
  rB1P2 = sqrt( (HSdiameter + ecc1*sen(theta))**2 + (ecc1*cos(theta))**2 )
  rB1B2 = HSdiameter
  rB1Q2 = sqrt( (HSdiameter - ecc*sen(theta))**2 + (ecc*cos(theta))**2 )
  rQ1P2 = sqrt( (HSdiameter - ecc2 - ecc1*sen(theta))**2 + (ecc1*cos(theta))**2 )
  rQ1B2 = HSdiameter - ecc2
  rQ1Q2 = sqrt( (HSdiameter - ecc2*(1 + sen(theta)))**2 + (ecc2*cos(theta))**2 )

  ws1s1 = computeOmega(patch1Radius, patch1Radius, rP1P2)
  ws1s2 = computeOmega(patch1Radius, patch2Radius, rP1Q2)
  ws2s2 = computeOmega(patch2Radius, patch2Radius, rQ1Q2)
  
  V = (
        eBB*    computeOmega(bigRadius,    bigRadius,    rB1B2) +
        eBs1*(  computeOmega(patch1Radius, bigRadius,    rP1B2) +
                computeOmega(bigRadius,    patch1Radius, rB1P2) ) +
        eBs2*(  computeOmega(bigRadius,    patch2Radius, rB1Q2) +
                computeOmega(patch2Radius, bigRadius,    rQ1B2) ) +
        es1s1*  computeOmega(patch1Radius, patch1Radius, rP1P2) +
        es1s2*( computeOmega(patch1Radius, patch2Radius, rP1Q2) +
                computeOmega(patch2Radius, patch1Radius, rQ1P2) ) +
        es2s2*  computeOmega(patch2Radius, patch2Radius, rQ1Q2)
      )
  return V, ws1s1, ws1s2, ws2s2

def computePotentials():
  thetas = [ i*2.*pi/100. for i in range(101) ]
  Vh = []
  Vv = []
  wBB = []
  wBs1 = []
  wBs2 = []
  ws1s1 = []
  ws1s2 = []
  ws2s2 = []
  for theta in thetas:
    a, b, c, d = computePotentialVerticalParticle(theta)
    e, f, g, h = computePotentialHorizontalParticle(theta)

    Vv.append(a)
    wBB.append(b)
    wBs1.append(c)
    wBs2.append(d)
    Vh.append(e)
    ws1s1.append(f)
    ws1s2.append(g)
    ws2s2.append(h)

  return Vv, Vh, wBB, wBS, wSS, thetas

fBB = computeOmega(bigRadius, bigRadius, HSdiameter)
fBS = computeOmega(bigRadius, patchRadius, HSdiameter - ecc)
fSS = computeOmega(patchRadius, patchRadius, HSdiameter - 2*ecc)

print("volumes at contact:")
print( fBB, fBS, fSS)

if usageMode == "yes":
  eBB /= fBB
  eBS /= fBS
  eSS /= fSS
  cBB = eBB / emin
  cBS = eBS / emin
  cSS = eSS / emin
  print("scaled OUTPUT:")
elif usageMode == "no":
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

Vv, Vh, wBB, wBS, wSS, theta = computePotentials()

outputWns = open("wBB.txt", 'w')
outputWss = open("wBBscaled.txt", 'w')
outputPot = open("potential.txt", 'w')
for Vvi, Vhi, wBBi, wBSi, wSSi, thetai in zip(Vv, Vh, wBB, wBS, wSS, theta):
  outputWns.write(str(thetai).ljust(24) + str(wBBi).ljust(24) +
                  str(wBSi).ljust(24) + str(wSSi).ljust(24) + "\n")
  outputWss.write(str(thetai).ljust(24) + str(wBBi/fBB).ljust(24) +
                  str(wBSi/fBS).ljust(24) + str(wSSi/fSS).ljust(24) + "\n")
  outputPot.write(str(thetai).ljust(24) + str(Vvi).ljust(24) + str(Vhi).ljust(24) +
                  str(Vvi/emin).ljust(24) + str(Vhi/emin).ljust(24) + "\n")
#for Vvi, Vhi, wBBi, wBSi, wSSi, thetai in zip(reversed(Vv), reversed(Vh), reversed(wBB), reversed(wBS), reversed(wSS), theta):
#  outputWns.write(str(thetai+0.5*pi).ljust(24) + str(wBBi).ljust(24) +
#                  str(wBSi).ljust(24) + str(wSSi).ljust(24) + "\n")
#  outputWss.write(str(thetai+0.5*pi).ljust(24) + str(wBBi/fBB).ljust(24) +
#                  str(wBSi/fBS).ljust(24) + str(wSSi/fSS).ljust(24) + "\n")
#  outputPot.write(str(thetai+0.5*pi).ljust(24) + str(Vvi).ljust(24) + str(Vhi).ljust(24) +
#                  str(Vvi/emin).ljust(24) + str(Vhi/emin).ljust(24) + "\n")
