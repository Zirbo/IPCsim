#! /usr/bin/python3

# writes a configuration with a plane structure of ipcs.

import argparse
from math import cos, sin, sqrt, pi
from numpy.random import ranf

parser = argparse.ArgumentParser(description='Creates a LAMMPS starting configuration.')
parser.add_argument('nIPCsSide', metavar='N', type=int, help='cubic root of the number of particles')
parser.add_argument('density', metavar='d', type=float, help='density')
parser.add_argument('ecc', metavar='e', type=float, help='eccentricity off the IPCs')
args = parser.parse_args()

outputFile = open('startingConfiguration.txt','w')
args.nIPCsPlane = args.nIPCsSide**2
args.nIPCs = args.nIPCsSide**3
args.side = (args.nIPCs/args.density)**(1./3.)

print(args)


outputFile.write("# 3D starting configuration for LAMMPS created with a script available at")
outputFile.write("# https://gitlab.com/catnat/ipc_brownian_motion")
outputFile.write("\n")
outputFile.write("\n" + str(3*args.nIPCs).rjust(10) + " atoms")
outputFile.write("\n" + str(2*args.nIPCs).rjust(10) + " bonds")
outputFile.write("\n" + str(0).rjust(10) + " angles")
outputFile.write("\n" + str(0).rjust(10) + " dihedrals")
outputFile.write("\n")

outputFile.write("\n" + str(2).rjust(10) + " atom types")
outputFile.write("\n" + str(1).rjust(10) + " bond types")
outputFile.write("\n" + str(0).rjust(10) + " angle types")
outputFile.write("\n" + str(0).rjust(10) + " dihedral types")
outputFile.write("\n")

outputFile.write("\n" + str(0.0).rjust(10) + '{:3.8f}'.format(args.side).rjust(16) + " xlo xhi")
outputFile.write("\n" + str(0.0).rjust(10) + '{:3.8f}'.format(args.side).rjust(16) + " ylo yhi")
outputFile.write("\n" + str(0.0).rjust(10) + '{:3.8f}'.format(args.side).rjust(16) + " zlo zhi")
outputFile.write("\n")
outputFile.write("\nMasses")
outputFile.write("\n" + str(1).rjust(10) + str(2.0).rjust(10))
outputFile.write("\n" + str(2).rjust(10) + str(0.5).rjust(10))
outputFile.write("\n")
outputFile.write("\nAtoms")
outputFile.write("\n")
for iz in range(0, args.nIPCsSide):
    for iy in range(0, args.nIPCsSide):
        for ix in range(0, args.nIPCsSide):
            # ipc center
            i = 3*(ix + iy*args.nIPCsSide + iz*args.nIPCsPlane)
            x = (0.5 + ix + 0.1*ranf())*args.side
            y = (0.5 + iy + 0.1*ranf())*args.side
            z = (0.5 + iz + 0.1*ranf())*args.side
            outputFile.write("\n" + str(i).rjust(10) + str(1).rjust(10) +
                 '{:3.8f}'.format(x).rjust(16) +
                 '{:3.8f}'.format(y).rjust(16) +
                 '{:3.8f}'.format(z).rjust(16))
            # create random unit vector
            px = ranf() 
            py = ranf()
            pz = ranf()
            mod = (px**2 + py**2 + pz**2)**(0.5)
            px /= mod
            py /= mod
            pz /= mod
            # patches
            outputFile.write("\n" + str(i+1).rjust(10) + str(2).rjust(10) +
                 '{:3.8f}'.format(x + args.ecc*px).rjust(16) +
                 '{:3.8f}'.format(y + args.ecc*py).rjust(16) +
                 '{:3.8f}'.format(z + args.ecc*pz).rjust(16))
            outputFile.write("\n" + str(i+2).rjust(10) + str(2).rjust(10) +
                 '{:3.8f}'.format(x - args.ecc*px).rjust(16) +
                 '{:3.8f}'.format(y - args.ecc*py).rjust(16) +
                 '{:3.8f}'.format(z - args.ecc*pz).rjust(16))

outputFile.write("\n")
