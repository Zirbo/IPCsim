#! /usr/bin/python3

# writes a configuration with a single manner plane of ipcs, and other ipcs above it.

from math import cos, sin, sqrt, pi, floor
from numpy.random import ranf as ran

outputFile = open('startingstate.xyz','w')

NPx = 14           # particles in a side
NPy = 12           # particles in a side
spacing = 1.2      # distance between fluid particles
L   = 15           # scaling (box side)
ecc = 0.16         # eccentricity


ecc /= L
cos30 = sqrt(3)*.5
NFx = int(L/spacing)
NFy = int(L/spacing)
NFz = int(L/spacing) - 1
spacing /= L

def absolutePBC(z):
    return z - floor(z)

tempN = 3*(NPx*NPy + NFx*NFy*NFz)
print('Estimate of the number of particles: ', tempN)
outputFile.write(str(tempN) + '\n0.0000\n')

vm = sqrt(.3)
def vel():
    return vm*(2.*ran()-1.)/L

alpha = .91*2*pi
beta = .51*2*pi
gamma = .78*2*pi
p = [ [ ecc*cos(alpha), ecc*sin(alpha), 0. ] ,
      [ ecc*cos(beta),  ecc*sin(beta),  0. ] ,
      [ ecc*cos(gamma), ecc*sin(gamma), 0. ] ,
      [ 0.           ,  0.           ,  1. ] ]
wafer, fluid = 0,0


# wafer layer
z = 0.25/L
for ix in range(NPx):
    x = (1.5 + .5 + 1.0000000000001*cos30*ix)/L
    for iy in range(NPy):
        wafer += 1
        # ipc center
        j = 0 if (iy + (int((ix + 1)/2))%2)%2==0 else 1
        y = ( 1.5 + ( (.5 + 1.0000000000001*iy) if ix%2==0 else (1.0000000000001*iy) ) )/L
        vx = vy = vz = 0.
        outputFile.write('W\t' + str(x).ljust(20) + '\t' + str(y).ljust(20) + '\t' + str(z).ljust(20) + '\t')
        outputFile.write(str(vx).ljust(20) + '\t' + str(vy).ljust(20) + '\t' + str(vz).ljust(20) + '\n')
        # first patch
        px = x + p[j][0];    px = absolutePBC(px)
        py = y + p[j][1];    py = absolutePBC(py)
        pz = z + p[j][2];    pz = absolutePBC(pz)
        outputFile.write('P\t' + str(px).ljust(20) + '\t' + str(py).ljust(20) + '\t' + str(pz).ljust(20) + '\t')
        outputFile.write(str(vx).ljust(20) + '\t' + str(vy).ljust(20) + '\t' + str(vz).ljust(20) + '\n')
        # second patch
        px = x - p[j][0];    px = absolutePBC(px)
        py = y - p[j][1];    py = absolutePBC(py)
        pz = z - p[j][2];    pz = absolutePBC(pz)
        outputFile.write('P\t' + str(px).ljust(20) + '\t' + str(py).ljust(20) + '\t' + str(pz).ljust(20) + '\t')
        outputFile.write(str(vx).ljust(20) + '\t' + str(vy).ljust(20) + '\t' + str(vz).ljust(20) + '\n')

# square lattice above the plane

for iz in range(1, NFz + 1):
    z = (.5 + iz)*spacing
    for ix in range(NFx):
        x = (.5 + ix)*spacing
        for iy in range(NFy):
            fluid += 1
            # ipc center
            y = (.5 + iy)*spacing
            vx = vel()
            vy = vel()
            vz = vel()
            outputFile.write('F\t' + str(x).ljust(20) + '\t' + str(y).ljust(20) + '\t' + str(z).ljust(20) + '\t')
            outputFile.write(str(vx).ljust(20) + '\t' + str(vy).ljust(20) + '\t' + str(vz).ljust(20) + '\n')
            # first patch
            px = x + p[j][0];    px = absolutePBC(px)
            py = y + p[j][1];    py = absolutePBC(py)
            pz = z + p[j][2];    pz = absolutePBC(pz)
            outputFile.write('P\t' + str(px).ljust(20) + '\t' + str(py).ljust(20) + '\t' + str(pz).ljust(20) + '\t')
            outputFile.write(str(vx).ljust(20) + '\t' + str(vy).ljust(20) + '\t' + str(vz).ljust(20) + '\n')
            # second patch
            px = x - p[j][0];    px = absolutePBC(px)
            py = y - p[j][1];    py = absolutePBC(py)
            pz = z - p[j][2];    pz = absolutePBC(pz)
            outputFile.write('P\t' + str(px).ljust(20) + '\t' + str(py).ljust(20) + '\t' + str(pz).ljust(20) + '\t')
            outputFile.write(str(vx).ljust(20) + '\t' + str(vy).ljust(20) + '\t' + str(vz).ljust(20) + '\n')

print('wafer = ', wafer, 'fluid = ', fluid)
print('tot = ', wafer + fluid)
print('Resulting density = ', (wafer + fluid)/L**3)
print('NUMBER OF POINT PARTICLES = ', 3*(wafer + fluid) )
print('HC radius in reduced units = ', 1/L)
outputFile.seek(0)
outputFile.write(str(3*(wafer + fluid)))
