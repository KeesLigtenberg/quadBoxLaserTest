#!/bin/bash

RUN=$1
if [[ -z $2 ]]; then 
root -l <<EOF
.X laserDataFitter/fitLaserData.h+("run${RUN}/run${RUN}_events.root", "run${RUN}/fitted.root", "run${RUN}/align.dat");
EOF
else 
root -l <<EOF
.L laserDataFitter/fitLaserData.h+
fitAndUpdateAlignment("run${RUN}/run${RUN}_events.root", "run${RUN}/fitted.root", "run${RUN}/align.dat", $2 );
EOF
fi
