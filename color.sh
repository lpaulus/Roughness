if [ $# -ne 4 ]; then
    echo "bash color.sh purefe1_t7 SmoothRadius CurvatureSmooth AverageRadius"
    exit 1
fi
name=$1
smooth=$2
curvature=$3
average=$4
dir="results/$name/${smooth}/${curvature}/${average}"
root=$(pwd)
if [ ! -d $dir ]; then
    echo "Directory $dir does not exists. Run 'bash process.sh' with the same argument to compute the roughness"
    exit 1
fi
cd $dir
if [ ! -f original_kmax.off ]; then
    $root/build/Color normalised.off original_kmax.txt
    mv output.off original_kmax.off
fi
if [ ! -f original_kmax_sqrt.off ]; then
    $root/build/Color normalised.off original_kmax.txt -s
    mv output.off original_kmax_sqrt.off
fi
if [ ! -f smooth_kmax.off ]; then
    $root/build/Color smooth.off smooth_kmax.txt
    mv output.off smooth_kmax.off
fi
if [ ! -f smooth_kmax_sqrt.off ]; then
    $root/build/Color smooth.off smooth_kmax.txt -s
    mv output.off smooth_kmax_sqrt.off
fi
$root/build/Color normalised.off original_kav.txt
mv output.off original_kav.off
$root/build/Color normalised.off original_kav.txt -s
mv output.off original_kav_sqrt.off
$root/build/Color smooth.off smooth_kav.txt
mv output.off smooth_kav.off
$root/build/Color smooth.off smooth_kav.txt -s
mv output.off smooth_kav_sqrt.off
$root/build/Color normalised.off roughness.txt
mv output.off normalised_roughness.off
$root/build/Color normalised.off roughness.txt -s
mv output.off normalised_roughness_sqrt.off
$root/build/Color smooth.off roughness.txt
mv output.off smooth_roughness.off
$root/build/Color smooth.off roughness.txt -s
mv output.off smooth_roughness_sqrt.off
