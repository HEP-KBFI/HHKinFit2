# update cmake when necessary i.e lxplus
if [ "${HOSTNAME:0:6}" = lxplus ]; then
    source /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/cmake/3.28.3-6dd20ba899a48e46b39e41512ff0af64/etc/profile.d/init.sh
fi
mkdir -p build
cd build
make clean
if [ "${HOSTNAME:0:6}" = lxplus ]; then
    echo "lxplus"
    echo "${USER}"
    cmake .. -DPYTHON_LIBRARY_DIR="." -DPYTHON_EXECUTABLE="/usr/bin/python3" -DLXPLUS=ON -DBIND11LOC="/afs/cern.ch/user/t/${USER}/.local/lib/python3.9/site-packages/pybind11/share/cmake/pybind11/"
else
    cmake .. -DPYTHON_LIBRARY_DIR="." -DPYTHON_EXECUTABLE="/usr/bin/python3"
fi
make -j4
cp hhkinfit2*.so ../python/ 
