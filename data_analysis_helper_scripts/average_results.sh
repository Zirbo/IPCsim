#! /bin/bash

#../file_averager.py averageNumberOfNeighbours.out 2
#../file_averager.py g_r 2
#../file_averager.py meanSquaredDisplacement.out 2
../file_averager.py autocorrelations.out 3
../file_averager.py averageClusterSizes.out 2
../file_averager.py averageNOP.out 2
../clusterSize_averager.py averageNumberOfNeighbours.out
../file_averager.py g_r.out 2
../file_averager.py meanSquaredDisplacement.out 2
