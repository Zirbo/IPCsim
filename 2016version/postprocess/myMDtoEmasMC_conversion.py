#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# Converts mio, supposed a startingstate.xyz, i.e. an output from my MD sims
# and converts it to the format Emanuela uses in her MC.

mio=open('SaraLaVoltaBuona.xyz')
suo=open('SaraLaVoltaBuona.dat','w')

# Put the correct L, and also fill N as a check that the files are right...
L = 21.544346900318832
Nparticles = 12000
zero = 0.0

cacca = mio.readline()
cacca = mio.readline()
del(cacca)
suo.write(str(Nparticles)+'\n')
suo.write(str(zero).rjust(20)+str(zero).rjust(20)+str(zero).rjust(20)+'\n')
suo.write(str(L).rjust(20)+str(zero).rjust(20)+str(zero).rjust(20)+'\n')
suo.write(str(zero).rjust(20)+str(L).rjust(20)+str(zero).rjust(20)+'\n')
suo.write(str(zero).rjust(20)+str(zero).rjust(20)+str(L).rjust(20)+'\n')
for line in mio:
  a = line.split()
  if a[0] == 'C':
    Nr = 1.0;   Np = 19
  elif a[0] == 'P':
    Nr = .56;   Np = 4
  else:
    print('Errore!!!\n')
    exit()
  x = float(a[1])*L
  y = float(a[2])*L
  z = float(a[3])*L
  suo.write( '{:16.8f}'.format(x).ljust(20)+'{:16.8f}'.format(y).ljust(20)
  +'{:16.8f}'.format(z).ljust(20)+'{:3.2f}'.format(Nr).ljust(20)+str(Np)+'\n' )
