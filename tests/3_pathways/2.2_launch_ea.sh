#!/bin/bash

DATENOW=$( date )
echo "Started at: $DATENOW"


for script in magma_scripts/*.sh; do
  while [ $( ps -f -u $USER | grep 'magma' | wc -l ) -ge 20 ]; do sleep 1; done
  ./$script > ${script}.log &
done
wait




DATENOW=$( date )
echo "Finished at: $DATENOW"

