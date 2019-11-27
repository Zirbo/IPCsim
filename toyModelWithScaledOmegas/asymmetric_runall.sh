#! /bin/bash

mkdir -p models-asym

echo 3645n
./asymmetric_computeAngularPlots.py no 0.2  0.28  0.22  1.06934   -24.0160  -16.7584 216.652   141.341   177.666 1
mv potential.txt models-asym/potential_3645n.txt
mv wBB.txt models-asym/3645n_wBB.txt
mv wBBscaled.txt models-asym/3645n_wBBscaled.txt
echo

echo 45n
./asymmetric_computeAngularPlots.py no 0.2  0.22 0.22 .245728 -3.11694  -3.11694  21.2298 21.2298 21.2298 0.142495
mv potential.txt models-asym/potential_45n.txt
mv wBB.txt models-asym/45n_wBB.txt
mv wBBscaled.txt models-asym/45n_wBBscaled.txt
echo
