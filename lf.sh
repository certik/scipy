git fetch ondrej
git submodule update --init
mkdir lfortran-build/
cd lfortran-build/
LIBRARY_PATH="/home/harshita/Desktop/lfortran/src/runtime/"
FC=lfortran cmake \
  -DCMAKE_Fortran_FLAGS=--verbose \
  -DLFORTRAN_RUNTIME_LIBRARY_PATH=$LIBRARY_PATH \
  ..
make install
cp /home/harshita/Desktop/lfortran/src/runtime/liblfortran_runtime.*  $CONDA_PREFIX/lib
cd ../
python dev.py build
python dev.py test -t scipy.odr -v