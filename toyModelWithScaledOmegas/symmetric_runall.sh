#! /bin/bash

mkdir -p models

echo isoW
./symmetric_computeAngularPlots.py yes 0.2 0.22 -1 0 0 1
mv potential.txt models/potential_isoW.txt
echo
echo cPat
./symmetric_computeAngularPlots.py yes 0.2 0.22 0 0 -1 1
mv potential.txt models/potential_cPat.txt
echo
echo minIPC
./symmetric_computeAngularPlots.py yes 0.2 0.22 0 -1 2 1
mv potential.txt models/potential_minIPC.txt
echo
echo 45n
./symmetric_computeAngularPlots.py no 0.2 0.22 .245728 -3.11694 21.2298 .142495
mv potential.txt models/potential_45n.txt
echo
echo 45n-0
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.1 -1.1 6.1 1.0
mv potential.txt models/potential_45n-0.txt
echo
echo 45n-1
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.0 -1.0 6.0 1.0
mv potential.txt models/potential_45n-1.txt
echo
echo 45n-2
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.1 -1.1 2.1 1.0
mv potential.txt models/potential_45n-2.txt
echo
echo 45c
./symmetric_computeAngularPlots.py no 0.2 0.22 3.18802 -24.3562 58.9717 1.02726
mv potential.txt models/potential_45c.txt
echo
echo 45c-0
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.2 -1.2 2.4 1.0
mv potential.txt models/potential_45c-0.txt
echo
echo 45c-I
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.0 -1.0 2.2 1.0
mv potential.txt models/potential_45c-1.txt
echo
echo 45c-II
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.2 -1.2 2.2 1.0
mv potential.txt models/potential_45c-2.txt

gnuplot symmetric_plotmodels.gpi
#gnuplot symmetric_plotmodelsPNG.gpi
mv *.eps models/
