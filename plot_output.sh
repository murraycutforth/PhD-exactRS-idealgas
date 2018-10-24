#!/bin/bash

for file in TTC1 TTC2 TTC3 TTC4 TTC5
do
	gnuplot -c epsplotsoln.gnu $file
done
