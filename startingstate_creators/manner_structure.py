#! /usr/bin/python3

# writes a configuration with a plane structure of ipcs.

from math import cos, sin, sqrt, pi
from numpy.random import ranf as ran

outputFile = open('startingstate.xyz','w')

Nx = 14            # particles in a side
Ny = 24            # particles in a side
Nz = 6             # # of planes
Lx = 12            # scaling (box side)
Ly = 24            # scaling (box side)
Lz = 12            # scaling (box side)
nIPCs = Nz*(Nx*Ny+int(Nx/2)*int((Ny+4)/2) )
ecc = 0.22         # eccentricity
cos30 = sqrt(3)*.5

outputFile.write(str(3*nIPCs) + '\n')
outputFile.write('Lx = {}    Ly = {}    Lz = {}\n'.format(Lx, Ly, Lz))

vm = sqrt(.3)
def vel():
#  return vm*(2.*ran()-1.)/L
  return 0.

alpha = .45*pi
beta = .93*pi
p = [ [ ecc*cos(alpha), ecc*sin(alpha), 0. ] ,
      [ ecc*cos(beta),  ecc*sin(beta),  0. ] ]
wafer, schokolade = 0,0


for iz in range(Nz):
  # wafer layer
  z = .5+2.0000000000001*iz
  for ix in range(Nx):
    x = .5+1.0000000000001*cos30*ix
    for iy in range(Ny):
      j = 0 if (iy+(int((ix+1)/2))%2)%2==0 else 1
      y = (.5+1.0000000000001*iy) if ix%2==0 else (1.0000000000001*iy)
      vx = vel()
      vy = vel()
      vz = vel()
      outputFile.write('W\t'+str(x).ljust(20)+'\t'+str(y).ljust(20)+'\t'+str(z).ljust(20)+'\t')
      outputFile.write(str(vx).ljust(20)+'\t'+str(vy).ljust(20)+'\t'+str(vz).ljust(20)+'\n')
      outputFile.write('P\t'+str(x+p[j][0]).ljust(20)+'\t'+str(y+p[j][1]).ljust(20)+
           '\t'+str(z+p[j][2]).ljust(20)+'\t')
      outputFile.write(str(vx).ljust(20)+'\t'+str(vy).ljust(20)+'\t'+str(vz).ljust(20)+'\n')
      outputFile.write('P\t'+str(x-p[j][0]).ljust(20)+'\t'+str(y-p[j][1]).ljust(20)+
           '\t'+str(z-p[j][2]).ljust(20)+'\t')
      outputFile.write(str(vx).ljust(20)+'\t'+str(vy).ljust(20)+'\t'+str(vz).ljust(20)+'\n')
      wafer += 1
for iz in range(Nz):
  # chocolate layer
  z = 1.5+2.0000000000001*iz
  for ix in range(0,Nx,2):
    x = .5+1.0000000000001*cos30*(ix+0.5)
    #x = (1.+1.0000000000001*(ix+0.5))/L
 #   for iy in range(0,Ny,2):
 #     y = (1.0000000000001*iy)/L
    for iy in range(0,Ny+4,2):     # golden
      y = 1.6 + 0.8*iy                    # golden
#    for iy in range(0,Ny,1):       # denser
#      y = iy/L                     # denser
#    for iy in range(0,Ny+1,3):       # less denser
#      y = (.8*iy)/L                  # less denser
      vx = vel()
      vy = vel()
      vz = vel()
      outputFile.write('S\t'+str(x).ljust(20)+'\t'+str(y).ljust(20)+'\t'+str(z).ljust(20)+'\t')
      outputFile.write(str(vx).ljust(20)+'\t'+str(vy).ljust(20)+'\t'+str(vz).ljust(20)+'\n')
      outputFile.write('Q\t'+str(x).ljust(20)+'\t'+str(y).ljust(20)+'\t'+str(z+ecc).ljust(20)+'\t')
      outputFile.write(str(vx).ljust(20)+'\t'+str(vy).ljust(20)+'\t'+str(vz).ljust(20)+'\n')
      outputFile.write('Q\t'+str(x).ljust(20)+'\t'+str(y).ljust(20)+'\t'+str(z-ecc).ljust(20)+'\t')
      outputFile.write(str(vx).ljust(20)+'\t'+str(vy).ljust(20)+'\t'+str(vz).ljust(20)+'\n')
      schokolade += 1
print('wafer =',wafer,'schokolade =',schokolade)
print('tot =',wafer+schokolade, 'ratio =',schokolade/(wafer+schokolade))
print('NUMBER OF POINT PARTICLES =',nIPCs,"=",wafer+schokolade )
print('density = ', nIPCs/(Lx*Ly*Lz))
#outputFile.seek(0)
#outputFile.write(str(3*(wafer+schokolade)))
