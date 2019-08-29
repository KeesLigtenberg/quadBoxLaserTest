#!/bin/bash

if [ ! -d run$1 ]; then
	mkdir run$1
fi

scp gastpx3@arawana:{/localstore/TPX3/laser_setup_software/run$1.root,/localstore2/TPX3/DATA/CHIP0/Run$1/stage_positions.txt} ./run$1/

cp align.dat run$1

