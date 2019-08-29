#!/bin/bash

RUN=$1
root -l <<EOF
.L eventBuilder/buildEvent.h+
convertToTree("run${RUN}/run${RUN}.root", "run${RUN}/stage_positions.txt", "run${RUN}/run${RUN}_events.root");
EOF
