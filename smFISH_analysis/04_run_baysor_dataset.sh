#! /usr/bin/bash
DIR=$(pwd)
for i in */
do
    echo Processing: $i
    DIR_INPUT="$DIR/$i"
    DIR_OUTPUT="$DIR/${i}output"
    mkdir -p $DIR_OUTPUT
    echo Output directory: $DIR_OUTPUT
    echo Initializing Baysor...
    baysor run -s 25 -o $DIR_OUTPUT ${DIR_INPUT}spots.csv ${DIR_INPUT}segmentation.tiff
done