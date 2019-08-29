#!/bin/bash

for i in {951..954}; do
echo "=========="
echo "RUN $i"
echo "=========="
echo " "
. doScp.sh $i
./buildEvent.sh $i 
./runFitted.sh $i 10 
done
