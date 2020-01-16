#! /usr/bin/python3

# writes a configuration with a plane structure of ipcs.

import argparse
from math import cos, sin, sqrt, pi, floor
from numpy.random import ranf

helpString = """Creates a LAMMPS starting configuration with a fully formed wafer layer structure.\n
Sample values for a cubic box:      14 12 6 12.4 12.4 12.0 0.22\n
Sample values for an elongated box: 14 24 6 12.4 24.8 12.0 0.22\n"""

parser = argparse.ArgumentParser(description=helpString)
parser.add_argument('particlePerSideX', metavar='nPx', type=int, help='number of IPCs in the X side')
parser.add_argument('particlePerSideY', metavar='nPy', type=int, help='number of IPCs in the Y side')
parser.add_argument('numberOfLayersZ', metavar='nLz', type=int, help='number of IPC layers in the Z side')
parser.add_argument('boxSideX', metavar='Lx', type=float, help='size of the simulation box side base (x)')
parser.add_argument('boxSideY', metavar='Ly', type=float, help='size of the simulation box side base (y)')
parser.add_argument('boxSideZ', metavar='Lz', type=float, help='height of the simulation box side (z)')
parser.add_argument('ecc', metavar='e', type=float, help='eccentricity of the IPCs')
args = parser.parse_args()
print(args)

outputFile = open('IPC_startingstate_manner.txt','w')

Lx = args.boxSideX
Ly = args.boxSideY
Lz = args.boxSideZ
ecc = args.ecc
nWaferX = args.particlePerSideX
nWaferY = args.particlePerSideY
nWaferZ = args.numberOfLayersZ
nFluidX = int(nWaferX/2)
nFluidY = int(nWaferY/1.6)
nFluidZ = nWaferZ

totalWaferIPCs = nWaferX * nWaferY * nWaferZ
totalChocoIPCs = nFluidX * nFluidY * nFluidZ
nIPCs = totalWaferIPCs + totalChocoIPCs

def absolutePBCx(x):
    return x - Lx*floor(x/Lx)

def absolutePBCy(y):
    return y - Ly*floor(y/Ly)

def absolutePBCz(z):
    return z - Lz*floor(z/Lz)

alpha = .45*pi
beta = .93*pi
cos30 = sqrt(3)*.5
p = [ [ ecc*cos(alpha), ecc*sin(alpha), 0. ] ,
      [ ecc*cos(beta),  ecc*sin(beta),  0. ] ]


outputFile.write("# 3D starting configuration for LAMMPS created with a script available at\n")
outputFile.write("# https://github.com/Zirbo/IPCsim/tree/master/lammps\n")
outputFile.write("# The wafer point particles are from 1 to " + str(3*totalWaferIPCs))

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
                        '{:3.8f}'.format(Lx).rjust(16) + "     xlo xhi")
outputFile.write("\n" + '{:3.8f}'.format(0.0).rjust(16) +
                        '{:3.8f}'.format(Ly).rjust(16) + "     ylo yhi")
outputFile.write("\n" + '{:3.8f}'.format(0.0).rjust(16) +
                        '{:3.8f}'.format(Lz).rjust(16) + "     zlo zhi")

outputFile.write("\n")
outputFile.write("\nMasses")
outputFile.write("\n#  atomtype, mass")
outputFile.write("\n" + str(1).rjust(10) + str(2.0).rjust(10))
outputFile.write("\n" + str(2).rjust(10) + str(0.5).rjust(10))

outputFile.write("\n")
outputFile.write("\nAtoms")
outputFile.write("\n#   atom-ID    mol-ID   atom-type    charge    x               y               z")


# wafer layer
waferIPCs = 0
for iz in range(nWaferZ):
    z = 0.5 + 2.0000000000001*iz
    for ix in range(nWaferX):
        x = 0.5 + 1.0000000000001*cos30*ix
        for iy in range(nWaferY):
            waferIPCs += 1
            atomNumber = (waferIPCs - 1)*3 + 1
            # ipc center
            j = 0 if (iy + (int((ix + 1)/2))%2)%2==0 else 1
            y = 0.5 + ( (.5 + 1.0000000000001*iy) if ix%2==0 else (1.0000000000001*iy) )
            outputFile.write("\n" + str(atomNumber).rjust(10) +
                  str(waferIPCs).rjust(10) +
                  str(1).rjust(10) +
                  str(-1.).rjust(10) +
                 '{:3.8f}'.format(x).rjust(16) +
                 '{:3.8f}'.format(y).rjust(16) +
                 '{:3.8f}'.format(z).rjust(16) )
            # first patch
            px = x + p[j][0];    px = absolutePBCx(px)
            py = y + p[j][1];    py = absolutePBCy(py)
            pz = z + p[j][2];    pz = absolutePBCz(pz)
            atomNumber += 1
            outputFile.write("\n" + str(atomNumber).rjust(10) +
                  str(waferIPCs).rjust(10) +
                  str(2).rjust(10) +
                  str(0.5).rjust(10) +
                 '{:3.8f}'.format(px).rjust(16) +
                 '{:3.8f}'.format(py).rjust(16) +
                 '{:3.8f}'.format(pz).rjust(16) )
            # second patch
            px = x - p[j][0];    px = absolutePBCx(px)
            py = y - p[j][1];    py = absolutePBCy(py)
            pz = z - p[j][2];    pz = absolutePBCz(pz)
            atomNumber += 1
            outputFile.write("\n" + str(atomNumber).rjust(10) +
                  str(waferIPCs).rjust(10) +
                  str(2).rjust(10) +
                  str(0.5).rjust(10) +
                 '{:3.8f}'.format(px).rjust(16) +
                 '{:3.8f}'.format(py).rjust(16) +
                 '{:3.8f}'.format(pz).rjust(16) )
# chocolate layer
chocoIPCs = 0
for iz in range(nFluidZ):
    z = 1.5 + 2.0000000000001*iz
    for ix in range(nFluidX):
        x = 0.6 + 1.0000000000001*cos30*(2*ix+0.5)
        for iy in range(nFluidY):
            chocoIPCs += 1
            atomNumber = 3*totalWaferIPCs + (chocoIPCs - 1)*3 + 1
            y = .6 + 1.6*iy
            outputFile.write("\n" + str(atomNumber).rjust(10) +
                  str(totalWaferIPCs + chocoIPCs).rjust(10) +
                  str(1).rjust(10) +
                  str(-1.).rjust(10) +
                 '{:3.8f}'.format(x).rjust(16) +
                 '{:3.8f}'.format(y).rjust(16) +
                 '{:3.8f}'.format(z).rjust(16) )
            # first patch
            px = x;              px = absolutePBCx(px)
            py = y;              py = absolutePBCy(py)
            pz = z + ecc;        pz = absolutePBCz(pz)
            atomNumber += 1
            outputFile.write("\n" + str(atomNumber).rjust(10) +
                  str(totalWaferIPCs + chocoIPCs).rjust(10) +
                  str(2).rjust(10) +
                  str(0.5).rjust(10) +
                 '{:3.8f}'.format(px).rjust(16) +
                 '{:3.8f}'.format(py).rjust(16) +
                 '{:3.8f}'.format(pz).rjust(16) )
            # second patch
            px = x;              px = absolutePBCx(px)
            py = y;              py = absolutePBCy(py)
            pz = z - ecc;        pz = absolutePBCz(pz)
            atomNumber += 1
            outputFile.write("\n" + str(atomNumber).rjust(10) +
                  str(totalWaferIPCs + chocoIPCs).rjust(10) +
                  str(2).rjust(10) +
                  str(0.5).rjust(10) +
                 '{:3.8f}'.format(px).rjust(16) +
                 '{:3.8f}'.format(py).rjust(16) +
                 '{:3.8f}'.format(pz).rjust(16) )

assert waferIPCs == totalWaferIPCs, '{} != {}'.format(waferIPCs, totalWaferIPCs)
assert chocoIPCs == totalChocoIPCs, '{} != {}'.format(chocoIPCs, totalChocoIPCs)
assert chocoIPCs + waferIPCs == nIPCs, '{} != {}'.format(chocoIPCs + waferIPCs, nIPCs)


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
