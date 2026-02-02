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
	echo "Please enter CONTROL $(( ${i}+1 )) sample: "
	read control
	controls+=(${control})
	echo
done


declare -a mutants
for (( i=0; i<$num_mutants; i++ ))
do
	echo "Please enter MUTANT $(( ${i}+1 )) sample: "
	read mutant
	mutants+=(${mutant})
	echo
done

echo "Please enter a summit bed file you would like to use for reference points: "
read tmp
bedFile=${tmp}

echo
echo "Please enter the histone modification of these samples: "
read modification

echo
directory=${modification}
mkdir ${directory}

echo ${#mutants[@]}

matrix_string=
heatmap_string=''

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

echo "computeMatrix scale-regions -p max -R ${bedFile} -S ${matrix_string} --upstream 500 --downstream 500 --regionBodyLength 2000 -o ${directory}matrix_${modification}.gz --missingDataAsZero --skipZeros"

computeMatrix scale-regions \
	-p max \
	-R ${bedFile} \
	-S ${matrix_string} \
	--upstream 500 \
	--downstream 500 \
	--regionBodyLength 2000 \
	-o ${directory}matrix_${modification}.gz \
	--missingDataAsZero --skipZeros


echo "plotHeatmap -m ${directory}matrix_${modification}.gz --colorMap 'Blues' --heatmapWidth 8 --heatmapHeight 64 --verbose --samplesLabel ${heatmap_string} -out ${directory}${modification}.png"

plotHeatmap -m ${directory}matrix_${modification}.gz --colorMap 'Blues' --heatmapWidth 8 --heatmapHeight 64 --verbose --samplesLabel ${heatmap_string} -out ${directory}${modification}.png
