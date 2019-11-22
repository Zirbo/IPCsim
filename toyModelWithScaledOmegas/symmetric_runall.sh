#! /bin/bash

mkdir -p models

echo 45n
./symmetric_computeAngularPlots.py no 0.2 0.22 .245728 -3.11694 21.2298 .142495
mv potential.txt models/potential_45n.txt
echo
echo 45c
./symmetric_computeAngularPlots.py no 0.2 0.22 3.18802 -24.3562 58.9717 1.02726
mv potential.txt models/potential_45c.txt
echo
echo 1_0_0
./symmetric_computeAngularPlots.py yes 0.2 0.22 -1 0 0 1
mv potential.txt models/potential_1_0_0.txt
echo
echo 0_0_1
./symmetric_computeAngularPlots.py yes 0.2 0.22 0 0 -1 1
mv potential.txt models/potential_0_0_1.txt
echo
echo 0_1_2
./symmetric_computeAngularPlots.py yes 0.2 0.22 0 -1 2 1
mv potential.txt models/potential_0_1_2.txt
echo
echo 01_11_61
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.1 -1.1 6.1 1.0
mv potential.txt models/potential_01_11_61.txt
echo
echo 02_12_24
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.2 -1.2 2.4 1.0
mv potential.txt models/potential_02_12_24.txt
