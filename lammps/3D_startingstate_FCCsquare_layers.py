#! /usr/bin/python3

# writes a configuration with a plane structure of ipcs.

import argparse
from math import cos, sin, sqrt, pi, floor
from numpy.random import ranf

helpString = """Creates a LAMMPS starting configuration with a single square plane.\n
Suggested values for a cubic box: 12 1.1 1.1 12.1 12.1 0.16\n
Suggested values for an elongates box for gravity experiments: 14 12 1.2 sz 12.4 Lz 0.22\n"""

parser = argparse.ArgumentParser(description=helpString)
parser.add_argument('particlePerSide', metavar='nPx', type=int, help='number of particles in the X-Y plane')
parser.add_argument('spacing', metavar='s', type=float, help='spacing of the fluid in the x-y plane')
parser.add_argument('spacingZ', metavar='sz', type=float, help='spacing of the fluid in the z direction')
parser.add_argument('boxSide', metavar='L', type=float, help='size of the simulation box side base (x-y)')
parser.add_argument('boxSideZ', metavar='Lz', type=float, help='height of the simulation box side (z)')
parser.add_argument('ecc', metavar='e', type=float, help='eccentricity of the IPCs')
args = parser.parse_args()
print(args)

outputFile = open('IPC_startingstate_FCCsquare.txt','w')

L = args.boxSide
Lz = args.boxSideZ
spacing = args.spacing
spacingZ = args.spacingZ
ecc = args.ecc
nPlaneX = args.particlePerSide
nFluidX = int(L/spacing)
nFluidY = int(L/spacing)
nFluidZ = int(Lz/spacingZ) - 1

print(nPlaneX, spacing, L, ecc, nFluidX, nFluidY, nFluidZ)

nIPCs = nPlaneX*nPlaneX + nFluidX*nFluidY*nFluidZ

def absolutePBC(x):
    return x - L*floor(x/L)

def absolutePBCz(z):
    return z - Lz*floor(z/Lz)

p = [ [ ecc, 0., 0. ] ,
      [ 0., ecc, 0. ] ]


outputFile.write("# 3D starting configuration for LAMMPS created with a script available at\n")
outputFile.write("# https://github.com/Zirbo/IPCsim/tree/master/lammps\n")
outputFile.write("# The plane particles are from 1 to " + str(nPlaneX*nPlaneX))

outputFile.write("\n")
outputFile.write("\n" + str(3*nIPCs).rjust(16) + " atoms")
outputFile.write("\n" + str(2*nIPCs).rjust(16) + " bonds")
outputFile.write("\n" + str(  nIPCs).rjust(16) + " angles")
outputFile.write("\n")

outputFile.write("\n" + str(2).rjust(16) + " atom types")
outputFile.write("\n" + str(1).rjust(16) + " bond types")
outputFile.write("\n" + str(1).rjust(16) + " angle types")
outputFile.write("\n")

outputFile.write("\n" + '{:3.8f}'.format(0.0).rjust(16) +
                        '{:3.8f}'.format(L).rjust(16) + "     xlo xhi")
outputFile.write("\n" + '{:3.8f}'.format(0.0).rjust(16) +
                        '{:3.8f}'.format(L).rjust(16) + "     ylo yhi")
outputFile.write("\n" + '{:3.8f}'.format(0.0).rjust(16) +
                        '{:3.8f}'.format(Lz).rjust(16) + "    zlo zhi")

outputFile.write("\n")
outputFile.write("\nMasses")
outputFile.write("\n#  atomtype, mass")
outputFile.write("\n" + str(1).rjust(10) + str(2.0).rjust(10))
outputFile.write("\n" + str(2).rjust(10) + str(0.5).rjust(10))

outputFile.write("\n")
outputFile.write("\nAtoms")
outputFile.write("\n#   atom-ID    mol-ID   atom-type    charge    x               y               z")

# wafer layer
waferParticles = 0
z = 0.25
for ix in range(nPlaneX):
    x = 0.5 + 1.002*ix
    for iy in range(nPlaneX):
        waferParticles += 1
        atomNumber = (waferParticles - 1)*3 + 1
        # ipc center
        j = 0 if (ix + iy) % 2 == 0 else 1
        y = 0.5 + 1.002*iy
        outputFile.write("\n" + str(atomNumber).rjust(10) +
              str(waferParticles).rjust(10) +
              str(1).rjust(10) +
              str(-1.).rjust(10) +
             '{:3.8f}'.format(x).rjust(16) +
             '{:3.8f}'.format(y).rjust(16) +
             '{:3.8f}'.format(z).rjust(16) )
        # first patch
        px = x + p[j][0];    px = absolutePBC(px)
        py = y + p[j][1];    py = absolutePBC(py)
        pz = z + p[j][2];    pz = absolutePBCz(pz)
        atomNumber += 1
        outputFile.write("\n" + str(atomNumber).rjust(10) +
              str(waferParticles).rjust(10) +
              str(2).rjust(10) +
              str(0.5).rjust(10) +
             '{:3.8f}'.format(px).rjust(16) +
             '{:3.8f}'.format(py).rjust(16) +
             '{:3.8f}'.format(pz).rjust(16) )
        # second patch
        px = x - p[j][0];    px = absolutePBC(px)
        py = y - p[j][1];    py = absolutePBC(py)
        pz = z - p[j][2];    pz = absolutePBCz(pz)
        atomNumber += 1
        outputFile.write("\n" + str(atomNumber).rjust(10) +
              str(waferParticles).rjust(10) +
              str(2).rjust(10) +
              str(0.5).rjust(10) +
             '{:3.8f}'.format(px).rjust(16) +
             '{:3.8f}'.format(py).rjust(16) +
             '{:3.8f}'.format(pz).rjust(16) )

# square lattice above the plane
fluidParticles = 0
for iz in range(1, nFluidZ + 1):
    z = (.5 + iz)*spacingZ
    for ix in range(nFluidX):
        x = (.5 + ix)*spacing
        for iy in range(nFluidY):
            fluidParticles += 1
            atomNumber = waferParticles*3 + (fluidParticles - 1)*3 + 1
            # ipc center
            y = (.5 + iy)*spacing
            outputFile.write("\n" + str(atomNumber).rjust(10) +
                  str(waferParticles + fluidParticles).rjust(10) +
                  str(1).rjust(10) +
                  str(-1.).rjust(10) +
                 '{:3.8f}'.format(x).rjust(16) +
                 '{:3.8f}'.format(y).rjust(16) +
                 '{:3.8f}'.format(z).rjust(16) )
            # first patch
            px = x + p[j][0];    px = absolutePBC(px)
            py = y + p[j][1];    py = absolutePBC(py)
            pz = z + p[j][2];    pz = absolutePBCz(pz)
            atomNumber += 1
            outputFile.write("\n" + str(atomNumber).rjust(10) +
                  str(waferParticles + fluidParticles).rjust(10) +
                  str(2).rjust(10) +
                  str(0.5).rjust(10) +
                 '{:3.8f}'.format(px).rjust(16) +
                 '{:3.8f}'.format(py).rjust(16) +
                 '{:3.8f}'.format(pz).rjust(16) )
            # second patch
            px = x - p[j][0];    px = absolutePBC(px)
            py = y - p[j][1];    py = absolutePBC(py)
            pz = z - p[j][2];    pz = absolutePBCz(pz)
            atomNumber += 1
            outputFile.write("\n" + str(atomNumber).rjust(10) +
                  str(waferParticles + fluidParticles).rjust(10) +
                  str(2).rjust(10) +
                  str(0.5).rjust(10) +
                 '{:3.8f}'.format(px).rjust(16) +
                 '{:3.8f}'.format(py).rjust(16) +
                 '{:3.8f}'.format(pz).rjust(16) )


outputFile.write("\n")
outputFile.write("\nBonds")
outputFile.write("\n#  ID bond-type atom-1 atom-2")
for i in range(nIPCs):
    IDcenter = 3*i + 1
    IDpatch1 = 3*i + 2
    IDpatch2 = 3*i + 3
    outputFile.write("\n" + str(2*i+1).rjust(10) + str(1).rjust(10) +
                            str(IDcenter).rjust(10) + str(IDpatch1).rjust(10) )
    outputFile.write("\n" + str(2*i+2).rjust(10) + str(1).rjust(10) +
                            str(IDcenter).rjust(10) + str(IDpatch2).rjust(10) )

outputFile.write("\n")
outputFile.write("\nAngles")
outputFile.write("\n#  ID    angle-type atom-1 atom-2 atom-3  (atom-2 is the center atom in angle)")
for i in range(nIPCs):
    IDcenter = 3*i + 1
    IDpatch1 = 3*i + 2
    IDpatch2 = 3*i + 3
    outputFile.write("\n" + str(i+1).rjust(10) + str(1).rjust(10) +
                            str(IDpatch1).rjust(10) + str(IDcenter).rjust(10) +
                            str(IDpatch2).rjust(10) )
outputFile.write("\n")
