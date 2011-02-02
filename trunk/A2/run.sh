#!/bin/csh
pwd
date
./apf -n 150 -t 400 -x 1 -y 1
./apf -n 151 -t 400 -x 1 -y 2
./apf -n 400 -t 200 -x 1 -y 4
./apf -n 299 -t 100 -x 2 -y 1
./apf -n 400 -t 100 -x 4 -y 1
date
