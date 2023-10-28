# iterate over all files of scipy/special/cdflib/ and

set -ex

# Compile the Fortran files

for src in $(find /Users/pranavchiku/repos/scipy/scipy/optimize/lbfgsb_src -name "*.f"); do
    lfortran -c $src -o ${src%.f}.o --fixed-form --implicit-typing --implicit-interface --generate-object-code --use-loop-variable-after-loop --all-mangling --mangle-underscore --bindc-mangling --legacy-array-sections
    cp ${src%.f}.o ./
done

lfortran -c ./scipy/optimize/lbfgsb_src/linpack.f -o linpack.o --fixed-form --implicit-typing --implicit-interface --generate-object-code --use-loop-variable-after-loop --all-mangling --mangle-underscore --bindc-mangling --legacy-array-sections --rtlib

# Link the object files into a shared library using clang

clang -shared -fPIC -o scipy/optimize/lbfgsb_src.so *.o -L/Users/pranavchiku/repos/lfortran/src/runtime -llfortran_runtime -llapack -lblas

mkdir -p $CONDA_PREFIX/lib/scipy/optimize/
cp scipy/optimize/lbfgsb_src.so $CONDA_PREFIX/lib/scipy/optimize/

# Copy the runtime library

# For Mac
cp ~/repos/lfortran/src/runtime/liblfortran_runtime.dylib $CONDA_PREFIX/lib/
cp ~/repos/lfortran/src/runtime/liblfortran_runtime.0.dylib $CONDA_PREFIX/lib/

# For Linux
# cp ~/repos/lfortran/src/runtime/liblfortran_runtime.so $CONDA_PREFIX/lib/

# Build

python dev.py build

# Setup the package

# For Mac
mkdir -p ./build-install/lib/python3.11/site-packages/scipy/optimize/
cp scipy/optimize/lbfgsb_src.so ./build-install/lib/python3.11/site-packages/scipy/optimize/


# For linux

#mkdir -p ./build-install/lib/python3.11/scipy/special/
#cp scipy/special/cdflib.so ./build-install/lib/python3.11/site-packages/scipy/special/

# Test

python dev.py test -t -s scipy.optimize -v
