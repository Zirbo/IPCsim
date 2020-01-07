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

def computePotentialVerticalParticle(theta):
  rP1P2 = sqrt( (HSdiameter + ecc*sen(theta))**2 + (ecc*(1 - cos(theta)))**2 )
  rP1B2 = sqrt( ecc*2 + HSdiameter**2 )
  rP1Q2 = sqrt( (HSdiameter - ecc*sen(theta))**2 + (ecc*(1 + cos(theta)))**2 )
  rB1P2 = sqrt( (HSdiameter + ecc*sen(theta))**2 + (ecc*cos(theta))**2 )
  rB1B2 = HSdiameter
  rB1Q2 = sqrt( (HSdiameter - ecc*sen(theta))**2 + (ecc*cos(theta))**2 )
  rQ1P2 = sqrt( (HSdiameter + ecc*sen(theta))**2 + (ecc*(1 + cos(theta)))**2 )
  rQ1B2 = sqrt( ecc*2 + HSdiameter**2 )
  rQ1Q2 = sqrt( (HSdiameter - ecc*sen(theta))**2 + (ecc*(1 - cos(theta)))**2 )

  wBB = computeOmega(bigRadius,   bigRadius,   rB1B2)
  wBS = computeOmega(bigRadius,   patchRadius, rB1Q2)
  
  V = (
        eBB*( computeOmega(bigRadius,   bigRadius,   rB1B2) ) +
        eBS*( computeOmega(patchRadius, bigRadius,   rP1B2) +
              computeOmega(bigRadius,   patchRadius, rB1P2) +
              computeOmega(bigRadius,   patchRadius, rB1Q2) +
              computeOmega(patchRadius, bigRadius,   rQ1B2) ) +
        eSS*( computeOmega(patchRadius, patchRadius, rP1P2) +
              computeOmega(patchRadius, patchRadius, rP1Q2) +
              computeOmega(patchRadius, patchRadius, rQ1P2) +
              computeOmega(patchRadius, patchRadius, rQ1Q2) )
      )
  return V, wBB, wBS

def computePotentialHorizontalParticle(theta):
  rP1P2 = sqrt( (HSdiameter + ecc*(1 + sen(theta)))**2 + (ecc*cos(theta))**2 )
  rP1B2 = ecc + HSdiameter
  rP1Q2 = sqrt( (HSdiameter + ecc*(1 - sen(theta)))**2 + (ecc*cos(theta))**2 )
  rB1P2 = sqrt( (HSdiameter + ecc*sen(theta))**2 + (ecc*cos(theta))**2 )
  rB1B2 = HSdiameter
  rB1Q2 = sqrt( (HSdiameter - ecc*sen(theta))**2 + (ecc*cos(theta))**2 )
  rQ1P2 = sqrt( (HSdiameter - ecc*(1 - sen(theta)))**2 + (ecc*cos(theta))**2 )
  rQ1B2 = HSdiameter - ecc
  rQ1Q2 = sqrt( (HSdiameter - ecc*(1 + sen(theta)))**2 + (ecc*cos(theta))**2 )

  wSS = computeOmega(patchRadius, patchRadius, rQ1Q2)
  
  V = (
        eBB*( computeOmega(bigRadius,   bigRadius,   rB1B2) ) +
        eBS*( computeOmega(patchRadius, bigRadius,   rP1B2) +
              computeOmega(bigRadius,   patchRadius, rB1P2) +
              computeOmega(bigRadius,   patchRadius, rB1Q2) +
              computeOmega(patchRadius, bigRadius,   rQ1B2) ) +
        eSS*( computeOmega(patchRadius, patchRadius, rP1P2) +
              computeOmega(patchRadius, patchRadius, rP1Q2) +
              computeOmega(patchRadius, patchRadius, rQ1P2) +
              computeOmega(patchRadius, patchRadius, rQ1Q2) )
      )
  return V, wSS

def computePotentials():
  thetas = [ i*pi/50. for i in range(51) ]
  Vh = []
  Vv = []
  wBB = []
  wBS = []
  wSS = []
  for theta in thetas:
    a, b, c = computePotentialVerticalParticle(theta)
    d, e = computePotentialHorizontalParticle(theta)

    Vh.append(d)
    Vv.append(a)
    wBB.append(b)
    wBS.append(c)
    wSS.append(e)

  return Vv, Vh, wBB, wBS, wSS, thetas

fBB = computeOmega(bigRadius, bigRadius, HSdiameter)
fBS = computeOmega(bigRadius, patchRadius, HSdiameter - ecc)
fSS = computeOmega(patchRadius, patchRadius, HSdiameter - 2*ecc)

#print("volumes at contact:")
#print( fBB, fBS, fSS)

if usageMode == "yes":
  eBB /= fBB
  eBS /= fBS
  eSS /= fSS
  cBB = eBB / emin
  cBS = eBS / emin
  cSS = eSS / emin
#  print("scaled OUTPUT:")
elif usageMode == "no":
  cBB = eBB * fBB / emin
  cBS = eBS * fBS / emin
  cSS = eSS * fSS / emin
#  print("multiplied OUTPUT:")
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
