#!/bin/bash

#read -p "WARNING: overwrites lowest energy structure [CTRL+C to cancel] " null
idx=${1}

touch tmp.sc
for files in $(ls part_?${idx}_*.sc)
do
    #tail -n+3 ${files} | sort -k2 -n | head -n 1 >> tmp.sc
    # get the top 10% instead
    #length=$(wc -l ${files} | awk '{print $1}')
    #slice_float=$(bc <<< "${length} * 0.1")
    #slice=${slice_float%.*}
    #tail -n+3 ${files} | sort -k2 -n | head -n ${slice} >> tmp.sc
    #cat ${files} >> tmp.sc
    tail -n+3 ${files} | sort -k2 -n >> tmp.sc
done

rm best_filenames${idx}.txt
touch best_filenames${idx}.txt
length=$(wc -l tmp.sc | awk '{print $1}')
slice_float=$(bc <<< "${length} * 0.1")
slice=${slice_float%.*}
tail -n+3 ${files} | sort -k2 -n | head -n ${slice} | awk '{print $(NF)}' >> best_filenames${idx}.txt
rm tmp.sc


