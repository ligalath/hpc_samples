#!/bin/bash
mpic++ oddevensort.cpp -o mpisort
mpiexec -n 4 ./mpisort