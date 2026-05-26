#!/bin/bash
target="$1"
name="$2"
folder="$3"
preprocessing_command="$4"
reso=8000

if [ "$target" = "empty" ]; then
    for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X
    do
        $preprocessing_command/preprocessing_mm10.sh $target $folder $reso chr$i $name &
    done
    wait
    exit
fi

# Batch 1: chr1-8
for i in 1 2 3 4 5 6 7 8; do
    $preprocessing_command/preprocessing_mm10.sh $target $folder $reso chr$i $name &
done
wait

# Batch 2: chr9-14
for i in 9 10 11 12 13 14; do
    $preprocessing_command/preprocessing_mm10.sh $target $folder $reso chr$i $name &
done
wait

# Batch 3: chr15-19, X
for i in 15 16 17 18 19 X; do
    $preprocessing_command/preprocessing_mm10.sh $target $folder $reso chr$i $name &
done
wait
