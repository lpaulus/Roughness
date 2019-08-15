if [ $# -ne 4 ]; then
    echo "bash process.sh purefe1_t7 SmoothRadius CurvatureSmooth AverageRadius"
    exit 1
fi
name=$1
smooth=$2
curvature=$3
average=$4
dir="results/$name/${smooth}_${curvature}_${average}"
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
$root/build/Roughness $root/data/$name.obj $smooth $curvature $average
if [ $? -ne 0 ]; then
    echo "Roughness failed"
    cd $root
    rm -r $dir
    exit $?
fi
cd $root
bash color.sh $name $smooth $curvature $average
