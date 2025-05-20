#!/bin/bash
SECONDS=0
input="data/halos_lcdm.txt"
output="result.hdf5"

ll=5.0

./run_mgs $ll $input $output

echo "Elapsed Time : $SECONDS sec"
exit 0
