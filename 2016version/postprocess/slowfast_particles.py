# -*- coding: utf-8 -*-

# see comments in slowfast_particles.sh

source = open("finalpos",'r')

data = (source.readline()).split()

farben = open("farben",'w')
fast = open("fastpart",'w')

less = 0 
more = 0
data.pop(0)
i=0

for x in data:
    i += 1
    ms = float(x)
    if ms <.5625:  # == (3/4)**2
        less += 1
        farben.write("B\n\n\n")
    else:
        more += 1
        farben.write("R\n\n\n")
        fast.write(str(i)+'\n')

print(less," mit  <x> < .75\n",more," mit <x> > .75\nsum = ", less+more)
