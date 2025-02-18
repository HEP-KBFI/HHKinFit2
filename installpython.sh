source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.24.08/x86_64-centos7-gcc48-opt/bin/thisroot.sh
mkdir -p build
cd build
make clean
cmake .. -DPYTHON_LIBRARY_DIR="." -DPYTHON_EXECUTABLE="/usr/bin/python3"
make -j4
