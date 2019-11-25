#! /usr/bin/python3

import argparse
from math import sqrt, pi, fabs, cos
from math import sin as sen

helpString = """Usage Modes: yes, no
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
parser.add_argument('es2s2', metavar='es2s2', type=float, help='')
parser.add_argument('es1s2', metavar='es1s2', type=float, help='')
parser.add_argument('emin', metavar='emin', type=float, help='')
##parser.add_argument('', metavar='', type=int, help='')
args = parser.parse_args()

usageMode = args.mode
ecc1 = args.ecc1
ecc2 = args.ecc2
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
  rQ1P2 = sqrt( (HSdiameter + ecc1*sen(theta))**2 + (ecc2 + ecc1*cos(theta))**2 )
  rQ1B2 = sqrt( ecc2*2 + HSdiameter**2 )
  rQ1Q2 = sqrt( (HSdiameter - ecc2*sen(theta))**2 + (ecc2*(1 - cos(theta)))**2 )

  wBB  = computeOmega(bigRadius, bigRadius,    rB1B2)
  wBs1 = computeOmega(bigRadius, patch1radius, rB1P2)
  wBs2 = computeOmega(bigRadius, patch2radius, rB1Q2)

  V = (
        eBB*    computeOmega(bigRadius,    bigRadius,    rB1B2) +
        eBs1*(  computeOmega(patch1radius, bigRadius,    rP1B2) +
                computeOmega(bigRadius,    patch1radius, rB1P2) ) +
        eBs2*(  computeOmega(bigRadius,    patch2radius, rB1Q2) +
                computeOmega(patch2radius, bigRadius,    rQ1B2) ) +
        es1s1*  computeOmega(patch1radius, patch1radius, rP1P2) +
        es1s2*( computeOmega(patch1radius, patch2radius, rP1Q2) +
                computeOmega(patch2radius, patch1radius, rQ1P2) ) +
        es2s2*  computeOmega(patch2radius, patch2radius, rQ1Q2)
      )
  return V, wBB, wBs1, wBs2

def computePotentialHorizontalParticle_P1left(theta):
  rP1P2 = sqrt( (HSdiameter + ecc1*(1 + sen(theta)))**2 + (ecc1*cos(theta))**2 )
  rP1B2 = ecc1 + HSdiameter
  rP1Q2 = sqrt( (HSdiameter + ecc1 - ecc2*sen(theta))**2 + (ecc2*cos(theta))**2 )
  rB1P2 = sqrt( (HSdiameter + ecc1*sen(theta))**2 + (ecc1*cos(theta))**2 )
  rB1B2 = HSdiameter
  rB1Q2 = sqrt( (HSdiameter - ecc2*sen(theta))**2 + (ecc2*cos(theta))**2 )
  rQ1P2 = sqrt( (HSdiameter - ecc2 + ecc1*sen(theta))**2 + (ecc1*cos(theta))**2 )
  rQ1B2 = HSdiameter - ecc2
  rQ1Q2 = sqrt( (HSdiameter - ecc2*(1 + sen(theta)))**2 + (ecc2*cos(theta))**2 )

  ws1s2 = computeOmega(patch2radius, patch1radius, rQ1P2)
  ws2s2 = computeOmega(patch2radius, patch2radius, rQ1Q2)

  V = (
        eBB*    computeOmega(bigRadius,    bigRadius,    rB1B2) +
        eBs1*(  computeOmega(patch1radius, bigRadius,    rP1B2) +
                computeOmega(bigRadius,    patch1radius, rB1P2) ) +
        eBs2*(  computeOmega(bigRadius,    patch2radius, rB1Q2) +
                computeOmega(patch2radius, bigRadius,    rQ1B2) ) +
        es1s1*  computeOmega(patch1radius, patch1radius, rP1P2) +
        es1s2*( computeOmega(patch1radius, patch2radius, rP1Q2) +
                computeOmega(patch2radius, patch1radius, rQ1P2) ) +
        es2s2*  computeOmega(patch2radius, patch2radius, rQ1Q2)
      )
  return V, ws1s2, ws2s2

def computePotentialHorizontalParticle_Q1left(theta):
  rP1P2 = sqrt( (HSdiameter + ecc1*(sen(theta) - 1))**2 + (ecc1*cos(theta))**2 )
  rP1B2 = HSdiameter - ecc1
  rP1Q2 = sqrt( (HSdiameter - ecc1 - ecc2*sen(theta))**2 + (ecc2*cos(theta))**2 )
  rB1P2 = sqrt( (HSdiameter + ecc1*sen(theta))**2 + (ecc1*cos(theta))**2 )
  rB1B2 = HSdiameter
  rB1Q2 = sqrt( (HSdiameter - ecc2*sen(theta))**2 + (ecc2*cos(theta))**2 )
  rQ1P2 = sqrt( (HSdiameter + ecc2 + ecc1*sen(theta))**2 + (ecc1*cos(theta))**2 )
  rQ1B2 = HSdiameter + ecc2
  rQ1Q2 = sqrt( (HSdiameter + ecc2*(1 - sen(theta)))**2 + (ecc2*cos(theta))**2 )

  ws1s1 = computeOmega(patch1radius, patch1radius, rP1P2)
  ws1s2 = computeOmega(patch1radius, patch2radius, rP1Q2)

  V = (
        eBB*    computeOmega(bigRadius,    bigRadius,    rB1B2) +
        eBs1*(  computeOmega(patch1radius, bigRadius,    rP1B2) +
                computeOmega(bigRadius,    patch1radius, rB1P2) ) +
        eBs2*(  computeOmega(bigRadius,    patch2radius, rB1Q2) +
                computeOmega(patch2radius, bigRadius,    rQ1B2) ) +
        es1s1*  computeOmega(patch1radius, patch1radius, rP1P2) +
        es1s2*( computeOmega(patch1radius, patch2radius, rP1Q2) +
                computeOmega(patch2radius, patch1radius, rQ1P2) ) +
        es2s2*  computeOmega(patch2radius, patch2radius, rQ1Q2)
      )
  return V, ws1s1, ws1s2

def computePotentials():
  thetas = [ i*2.*pi/100. for i in range(101) ]
  Vv = []
  VhP = []
  VhQ = []
  wBB = []
  wBs1 = []
  wBs2 = []
  ws1s1 = []
  ws1s2 = []
  ws2s2 = []
  for theta in thetas:
    a, b, c, d = computePotentialVerticalParticle(theta)
    e, f, g    = computePotentialHorizontalParticle_P1left(theta)
    h, i, j    = computePotentialHorizontalParticle_Q1left(theta)

    Vv.append(a)
    wBB.append(b)
    wBs1.append(c)
    wBs2.append(d)
    VhP.append(e)
    ws2s2.append(g)
    VhQ.append(h)
    ws1s1.append(i)
    ws1s2.append(f+j)

  return Vv, VhP, VhQ, wBB, wBs1, wBs2, ws1s1, ws1s2, ws2s2, thetas

fBB   = computeOmega(bigRadius,    bigRadius,    HSdiameter)
fBs1  = computeOmega(bigRadius,    patch1radius, HSdiameter - ecc1)
fBs2  = computeOmega(bigRadius,    patch2radius, HSdiameter - ecc2)
fs1s1 = computeOmega(patch1radius, patch1radius, HSdiameter - ecc1 - ecc1)
fs1s2 = computeOmega(patch1radius, patch2radius, HSdiameter - ecc1 - ecc2)
fs2s2 = computeOmega(patch2radius, patch2radius, HSdiameter - ecc2 - ecc2)

print("volumes at contact:")
print( fBB, fBs1, fBs2, fs1s1, fs1s2, fs2s2)

if usageMode == "yes":
  eBB   /= fBB
  eBs1  /= fBs1
  eBs2  /= fBs2
  es1s1 /= fs1s1
  es1s2 /= fs1s2
  es2s2 /= fs2s2
  cBB    = eBB   / emin
  cBs1   = eBs1  / emin
  cBs2   = eBs2  / emin
  cs1s1  = es1s1 / emin
  cs1s2  = es1s2 / emin
  cs2s2  = es2s2 / emin
  print("scaled OUTPUT:")
elif usageMode == "no":
  cBB   = eBB   * fBB   / emin
  cBs1  = eBs1  * fBs1  / emin
  cBs2  = eBs2  * fBs2  / emin
  cs1s1 = es1s1 * fs1s1 / emin
  cs1s2 = es1s2 * fs1s2 / emin
  cs2s2 = es2s2 * fs2s2 / emin
  print("multiplied OUTPUT:")
else:
  print("mode does not exist")
  exit()


print(str(cBB)[0:7].ljust(10)   + str(cBs1)[0:8].ljust(10)  + str(cBs2)[0:8].ljust(10))
print(str(cs1s1)[0:7].ljust(10) + str(cs1s2)[0:7].ljust(10) + str(cs2s2)[0:7].ljust(10))

Vv, VhP, VhQ, wBB, wBs1, wBs2, ws1s1, ws1s2, ws2s2, theta = computePotentials()

outputWns = open("wBB.txt", 'w')
outputWss = open("wBBscaled.txt", 'w')
outputPot = open("potential.txt", 'w')
for Vvi, VhPi, VhQi, wBBi, wBs1i, wBs2i, ws1s1i, ws1s2i, ws2s2i, thetai in zip(Vv, VhP, VhQ, wBB, wBs1, wBs2, ws1s1, ws1s2, ws2s2, theta):
#  outputWns.write(str(thetai).ljust(24) + str(wBBi).ljust(24) +
#                  str(wBs1i).ljust(24)  + str(wBs2i).ljust(24) +
#                  str(ws1s1i).ljust(24) + str(ws1s2i).ljust(24) +
#                  str(ws2s2i).ljust(24) + "\n")
#  outputWss.write(str(thetai).ljust(24)       + str(wBBi/fBB).ljust(24) +
#                  str(wBs1i/fBs1).ljust(24)   + str(wBs2i/fBs2).ljust(24) +
#                  str(ws1s1i/fs1s1).ljust(24) + str(ws1s2i/fs1s2).ljust(24) +
#                  str(ws2s2i/fs2s2).ljust(24) + "\n")
#  outputPot.write(str(thetai).ljust(24) +
#                  str(Vvi).ljust(24)       + str(VhPi).ljust(24) + str(VhQi).ljust(24) +
#                  str(Vvi/emin).ljust(24)  + str(VhPi/emin).ljust(24) +
#                  str(VhPi/emin).ljust(24) + "\n")
  outputWns.write(str(thetai).ljust(24) + str(wBBi).ljust(24) +
                  str(wBs1i + wBs2i).ljust(24) +
                  str(ws1s1i + ws2s2i).ljust(24) +
                  str(ws1s2i).ljust(24) + "\n")
  outputWss.write(str(thetai).ljust(24)       + str(wBBi/fBB).ljust(24) +
                  str(wBs1i/fBs1 + wBs2i/fBs2).ljust(24) +
                  str(ws1s1i/fs1s1 + ws2s2i/fs2s2).ljust(24) +
                  str(ws1s2i/fs1s2).ljust(24) + "\n")
  outputPot.write(str(thetai).ljust(24) + str(Vvi).ljust(24) +
                  str(VhPi).ljust(24) + str(VhQi).ljust(24) + "\n")
