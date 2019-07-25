#! /usr/bin/python3

# writes a configuration with a single manner plane of ipcs, and other ipcs above it.

from math import cos, sin, sqrt, pi
from numpy.random import ranf as ran

outputFile = open('startingstate.xyz','w')

NPx = 14           # particles in a side
NPy = 12           # particles in a side
spacing = 1.5      # distance between fluid particles
L   = 12           # scaling (box side)
ecc = 0.22         # eccentricity


ecc /= L
cos30 = sqrt(3)*.5
NFx = int(L/spacing)
NFy = int(L/spacing)
NFz = int(L/spacing) - 1
spacing /= L

tempN = 3*(NPx*NPy + NFx*NFy*NFz)
print('Estimate of the number of particles: ', tempN)
outputFile.write(str(tempN) + '\n0.0000\n')

vm = sqrt(.3)
def vel():
    return vm*(2.*ran()-1.)/L

alpha = .45*pi
beta = .93*pi
p = [ [ ecc*cos(alpha), ecc*sin(alpha), 0. ] ,
      [ ecc*cos(beta),  ecc*sin(beta),  0. ] ]
wafer, fluid = 0,0


# wafer layer
z = 0.
for ix in range(NPx):
    x = (.5 + 1.0000000000001*cos30*ix)/L
    for iy in range(NPy):
        j = 0 if (iy + (int((ix + 1)/2))%2)%2==0 else 1
        y = ((.5 + 1.0000000000001*iy) if ix%2==0 else (1.0000000000001*iy))/L
        vx = vy = vz = 0.
        # print the coordinates to the file
        outputFile.write('W\t' + str(x).ljust(20) + '\t' + str(y).ljust(20) + '\t' + str(z).ljust(20) + '\t')
        outputFile.write(str(vx).ljust(20) + '\t' + str(vy).ljust(20) + '\t' + str(vz).ljust(20) + '\n')
        outputFile.write('P\t' + str(x + p[j][0]).ljust(20) + '\t' + str(y + p[j][1]).ljust(20) + '\t' + str(z + p[j][2]).ljust(20) + '\t')
        outputFile.write(str(vx).ljust(20) + '\t' + str(vy).ljust(20) + '\t' + str(vz).ljust(20) + '\n')
        outputFile.write('P\t' + str(x-p[j][0]).ljust(20) + '\t' + str(y-p[j][1]).ljust(20) + '\t' + str(z-p[j][2]).ljust(20) + '\t')
        outputFile.write(str(vx).ljust(20) + '\t' + str(vy).ljust(20) + '\t' + str(vz).ljust(20) + '\n')
        # update the coutputFilenter
        wafer += 1

# square lattice above the plane

for iz in range(1, NFz + 1):
    z = (iz)*spacing
    for ix in range(NFx):
        x = (.5 + ix)*spacing
        for iy in range(NFy):
            y = (.5 + iy)*spacing
            vx = vel()
            vy = vel()
            vz = vel()
            # print the coordinates to the file
            outputFile.write('F\t' + str(x).ljust(20) + '\t' + str(y).ljust(20) + '\t' + str(z).ljust(20) + '\t')
            outputFile.write(str(vx).ljust(20) + '\t' + str(vy).ljust(20) + '\t' + str(vz).ljust(20) + '\n')
            outputFile.write('P\t' + str(x + p[j][0]).ljust(20) + '\t' + str(y + p[j][1]).ljust(20) + '\t' + str(z + p[j][2]).ljust(20) + '\t')
            outputFile.write(str(vx).ljust(20) + '\t' + str(vy).ljust(20) + '\t' + str(vz).ljust(20) + '\n')
            outputFile.write('P\t' + str(x-p[j][0]).ljust(20) + '\t' + str(y-p[j][1]).ljust(20) + '\t' + str(z-p[j][2]).ljust(20) + '\t')
            outputFile.write(str(vx).ljust(20) + '\t' + str(vy).ljust(20) + '\t' + str(vz).ljust(20) + '\n')
            # update the coutputFilenter
            fluid += 1

print('wafer =', wafer, 'fluid =', fluid)
print('tot =', wafer + fluid)
print('NUMBER OF POINT PARTICLES =', 3*(wafer + fluid) )
outputFile.seek(0)
outputFile.write(str(3*(wafer + fluid)))
