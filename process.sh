#!/bin/bash

input_sync () {
    if [ -f "$1" ]; then
        ln -s "$1"
    fi
}
output_sync () {
    if [ ! -f "$1/$2" ]; then
        mv "$2" "$1/$2"
        input_sync "$1/$2"
    fi
}

if [ $# -ne 4 ]; then
    echo "bash process.sh purefe1_t7 SmoothRadius CurvatureSmooth AverageRadius"
    exit 1
fi
name=$1
smooth=$2
curvature=$3
average=$4
dir="results/$name/${smooth}/${curvature}/${average}"
root=$(pwd)
if [ ! -f data/$name.obj ]; then
    echo "No $name.obj in data directory"
    exit 1
fi
if [ -d $dir ]; then
    echo "Directory $dir already exists. Run 'rm -r $dir' to get rid of the results."
    exit 1
fi
mkdir -p $dir
cd $dir
smooth_radius=$smooth
if [ -f "../../smooth.off" ]; then
    smooth_radius="-$smooth"
fi
input_sync "../../smooth.off"
$root/build/Roughness $root/data/$name.obj $smooth_radius $curvature $average
if [ $? -ne 0 ]; then
    echo "Roughness failed"
    cd $root
    rm -r $dir
    exit $?
fi
output_sync "../.." "smooth.off"
output_sync ".." "original_kmax.txt"
output_sync ".." "smooth_kmax.txt"
cd $root
input_sync "../original_kmax.off"
input_sync "../smooth_kmax.off"
input_sync "../original_kmax_sqrt.off"
input_sync "../smooth_kmax_sqrt.off"
bash color.sh $name $smooth $curvature $average
cd $dir
output_sync ".." "original_kmax.off"
output_sync ".." "smooth_kmax.off"
output_sync ".." "original_kmax_sqrt.off"
output_sync ".." "smooth_kmax_sqrt.off"
cd $root
