#! /usr/bin/python3

# writes a configuration with a single square plane of ipcs, and other ipcs above it.

from math import sqrt, floor
from numpy.random import ranf as ran

outputFile = open('startingstate.xyz','w')

NPx = 12           # particles in a side
NPy = 12           # particles in a side
spacing = 1.1      # distance between fluid particles
L   = 12.100       # scaling (box side)
ecc = 0.16         # eccentricity


ecc /= L
NFx = int(L/spacing)
NFy = int(L/spacing)
NFz = int(L/spacing) - 1
spacing /= L

def absolutePBC(z):
    return z - floor(z)

tempN = 3*(NPx*NPy + NFx*NFy*NFz)
print('Estimate of the number of particles: ', tempN)
outputFile.write(str(tempN) + '\n' + str(L) + '\t' + '0.0000\n')

vm = sqrt(.3)
def vel():
    return vm*(2.*ran()-1.)/L

p = [ [ ecc, 0., 0. ] ,
      [ 0., ecc, 0. ] ]
wafer, fluid = 0,0


# wafer layer
z = 0.25/L
for ix in range(NPx):
    x = ( .5 + 1.002*ix)/L
    for iy in range(NPy):
        wafer += 1
        # ipc center
        j = 0 if ( (ix + iy) % 2 == 0 ) else 1
        y = ( .5 + 1.002*iy )/L
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
