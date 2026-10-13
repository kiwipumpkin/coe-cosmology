#!/bin/bash
cd /workspace/clump_scan/sequential
P="python3 -u scenario_a.py outer"
($P random 1.5 10 101 > log_a_random_101.txt 2>&1) &
($P random 1.5 10 102 > log_a_random_102.txt 2>&1) &
($P random 1.5 10 103 > log_a_random_103.txt 2>&1) &
($P foam 1.5 15 201 > log_a_foam_201.txt 2>&1) &
($P foam 1.5 15 202 > log_a_foam_202.txt 2>&1) &
($P lattice 1.5 30 301 > log_a_lattice_301.txt 2>&1; $P random 1.5 6 601 octa > log_a_randomocta_601.txt 2>&1) &
($P pairs 0 200 401 > log_a_pairs_401.txt 2>&1) &
($P random 0.75 4 501 > log_a_random075_501.txt 2>&1) &
wait
echo ALLDONE
