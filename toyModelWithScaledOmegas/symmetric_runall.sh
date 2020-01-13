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
echo target
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.1 -1.1 6.1 1.0
mv potential.txt models/potential_target.txt
echo
echo PP1
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.0 -1.0 2.2 1.0
mv potential.txt models/potential_PP1.txt
echo
echo PP2
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.0 -1.0 4.0 1.0
mv potential.txt models/potential_PP2.txt
echo
echo PP3
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.0 -1.0 6.0 1.0
mv potential.txt models/potential_PP3.txt
echo
echo EE1
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.1 -1.1 2.1 1.0
mv potential.txt models/potential_EE1.txt
echo
echo EE2
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.2 -1.2 2.2 1.0
mv potential.txt models/potential_EE2.txt
echo
echo EE3
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.3 -1.3 2.3 1.0
mv potential.txt models/potential_EE3.txt
echo
echo EE4
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.4 -1.4 2.4 1.0
mv potential.txt models/potential_EE4.txt
echo
echo EE5
./symmetric_computeAngularPlots.py yes 0.2 0.22 0.5 -1.5 2.5 1.0
mv potential.txt models/potential_EE5.txt

gnuplot symmetric_plotmodels.gpi
for file in *.eps; do
  epspdf $file 
done
rm *.eps
mv *.pdf models/
