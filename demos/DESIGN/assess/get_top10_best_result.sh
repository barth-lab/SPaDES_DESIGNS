#!/bin/bash

length=$(ls *.pdb | wc -l | awk '{print $1}')
slice_float=$(bc <<< "${length} * 0.1")
slice=${slice_float%.*}

tail -n+3 score.sc | sort -k2 -n | head -n ${slice} >> top10.sc
