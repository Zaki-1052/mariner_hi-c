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
	echo "Please enter CONTROL $(( ${i}+1 )) sample (relative to /data/rs_256/genomewide_plots: "
	read control
	controls+=(/data/rs_256/genomewide_plots/${control})
	echo
done


declare -a mutants
for (( i=0; i<$num_mutants; i++ ))
do
	echo "Please enter MUTANT $(( ${i}+1 )) sample (relative to "/data/rs_256/genomewide_plots/"): "
	read mutant
	mutants+=(/data/rs_256/genomewide_plots/${mutant})
	echo
done

echo "Please enter a summit bed file you would like to use for reference points (relative to "/data/rs_256/genomewide_plots/"): "
read tmp
bedFile=/data/rs_256/genomewide_plots/${tmp}

echo
echo "Please enter the histone modification of these samples: "
read modification

echo
directory=/data/rs_256/genomewide_plots/${modification}/
mkdir ${directory}

echo ${#mutants[@]}

matrix_string=
heatmap_string=

if [[ ${num_controls} -gt ${num_mutants} ]]
then
	for i in "${!controls[@]}"
	do
		matrix_string="${matrix_string} ${controls[$i]}"
		heatmap_string="${heatmap_string} ${modification}_ctrl"
		if [[ i -lt ${#mutants[@]} ]]
		then
			matrix_string="${matrix_string} ${mutants[$i]}"
			heatmap_string="${heatmap_string} ${modification}_mut"
		fi
	done
else
	
	for i in "${!mutants[@]}"
	do
		if [[ i -lt ${#controls[@]} ]]
		then
			matrix_string="${matrix_string} ${controls[$i]}"
			heatmap_string="${heatmap_string} ${modification}_ctrl"
		fi
		matrix_string="${matrix_string} ${mutants[$i]}"
		heatmap_string="${heatmap_string} ${modification}_mut"
	done
fi

echo "creating matrix: computeMatrix reference-point --referencePoint center -R ${bedFile} -S${matrix_string} -a 1500 -b 1500 -o ${directory}matrix_${modification}.gz --missingDataAsZero"

computeMatrix reference-point -p max --referencePoint center -R ${bedFile} -S ${matrix_string} -a 5000 -b 5000 -o ${directory}matrix_${modification}.gz --missingDataAsZero

echo "plotting heatmap: plotHeatmap -m ${directory}matrix_${modification}.gz --samplesLabel ${heatmap_string} -out ${directory}${modification}.png"

plotHeatmap -m ${directory}matrix_${modification}.gz --colorMap 'Blues' --heatmapWidth 8 --heatmapHeight 64 --samplesLabel ${heatmap_string} -out ${directory}${modification}.png
