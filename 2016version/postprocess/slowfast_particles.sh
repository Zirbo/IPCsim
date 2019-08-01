#!/bin/bash

# run the program with source single_particle_msd.cpp. The output will
# contain the single particle mean squared displacement. The final
# positions will be sent to slowfast_particles.py that will produce:
# fastpart: the numbers of the fast particles, input for discerned_acs.cpp
# farben, used to modify the .xyz file to print in colors

tail -n 1 msd_singleparticles.anaout > finalpos
python3 slowfast_particles.py
rm finalpos
