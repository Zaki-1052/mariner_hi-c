#!/bin/bash

echo "Please enter the number of control samples"
read num_controls
echo

echo "Please enter the number of mutant samples"
read num_mutants
echo

declare -a controls
for (( i=0; i<$num_controls; i++ ))
do
    echo "Please enter CONTROL sample (FULL PATH): "
    read control
    controls+=("${control}")
    echo
done

declare -a mutants
for (( i=0; i<$num_mutants; i++ ))
do
    echo "Please enter MUTANT sample (FULL PATH): "
    read mutant
    mutants+=("${mutant}")
    echo
done

echo "Please enter a summit bed file you would like to use for reference points (FULL PATH): "
read bedFile

echo
echo "Please enter the histone modification of these samples: "
read modification

echo
echo "Please enter the directory you would like the heatmap in: "
read directory

matrix_string=""
heatmap_labels=""

# Interleave control and mutant samples for proper labeling
for (( i=0; i<${#controls[@]}; i++ ))
do
    matrix_string+=" ${controls[i]}"
    heatmap_labels+=" ${modification}_ctrl"
    if [[ i -lt ${#mutants[@]} ]]; then
        matrix_string+=" ${mutants[i]}"
        heatmap_labels+=" ${modification}_mut"
    fi
done

# Add any remaining mutant samples if more mutants than controls
for (( i=${#controls[@]}; i<${#mutants[@]}; i++ ))
do
    matrix_string+=" ${mutants[i]}"
    heatmap_labels+=" ${modification}_mut"
done

echo "Running computeMatrix..."
computeMatrix reference-point \
    --referencePoint center \
    -R "${bedFile}" \
    -S ${matrix_string} \
    -a 1500 -b 1500 \
    -o "${directory}matrix_${modification}.gz"

echo "Generating heatmap..."
plotHeatmap -m "${directory}matrix_${modification}.gz" \
    --samplesLabel ${heatmap_labels} \
    -out "${directory}${modification}.png"

echo "Done! Heatmap saved to ${directory}${modification}.png"

