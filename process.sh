if [ $# -ne 2 ]; then
    echo "Need 2 arguments"
    exit 1
fi
name=$1
eps=$2
dir="results/$name/$eps"
if [ ! -f data/$name.obj ]; then
    echo "No $name.obj in data directory"
    exit 1
fi
if [ -d $dir ]; then
    echo "Directory $dir already exists"
    exit 1
fi
mkdir -p $dir
cd $dir
../../../build/Roughness ../../../data/$name.obj $eps
if [ $? -ne 0 ]; then
    echo "Roughness failed"
    cd ../../../
    rm -r $dir
    exit $?
fi
../../../build/Color normalised.off original_kmax.txt
mv output.off original_kmax.off
../../../build/Color smooth.off smooth_kmax.txt
mv output.off smooth_kmax.off
../../../build/Color normalised.off original_kav.txt
mv output.off original_kav.off
../../../build/Color smooth.off smooth_kav.txt
mv output.off smooth_kav.off
../../../build/Color normalised.off roughness.txt
mv output.off normalised_roughness.off
../../../build/Color smooth.off roughness.txt
mv output.off smooth_roughness.off
