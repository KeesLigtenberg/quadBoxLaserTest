#!/bin/bash

RUN=$1
root -l <<EOF
.X laserDataFitter/fitLaserData.h+("run${RUN}/run${RUN}_events.root", "run${RUN}/run${RUN}_fitted.root", "run${RUN}/align.dat");
EOF
