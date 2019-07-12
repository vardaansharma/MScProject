#!/bin/bash

echo Running Geodesic algo
#./run -file ../datasets/PPI_Network.csv -skip 1 -a 1 -cols 2

echo Running Random edge
./run -file ../datasets/PPI_Network.csv -skip 1 -a 2 -cols 2 -quite

echo Running Spectral
./run -file ../datasets/PPI_Network.csv -skip 1 -a 3 -cols 2
