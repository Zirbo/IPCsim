#! /bin/bash

mkdir -p models-asym

echo 3645n
./asymmetric_computeAngularPlots.py no 0.2  0.28  0.22  1.06934   -24.0160  -16.7584 216.652   141.341   177.666 1
mv potential.txt models-asym/potential_3645n.txt
mv w* models-asym/
echo
