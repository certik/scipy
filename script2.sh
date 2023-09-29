# iterate over all files of scipy/special/amos/ and

set -ex

# Compile the Fortran files
for src in $(find /home/harshita/scipy/scipy/special/mach -name "*.f"); do
    gfortran -c $src -o ${src%.f}.o
    cp ${src%.f}.o ./
done

for src in $(find /home/harshita/scipy/scipy/special/amos -name "*.f"); do
    gfortran -c $src -o ${src%.f}.o 
    cp ${src%.f}.o ./
done

# lfortran -c ./scipy/special/amos/zabs.f -o zabs.o --fixed-form --implicit-typing --implicit-interface --generate-object-code --use-loop-variable-after-loop --all-mangling --mangle-underscore --bindc-mangling --rtlib

# Link the object files into a shared library using clang



gfortran -c ./scipy/special/amos/zbesh.f -o zbesh.o

# Link the object files into a shared library using clang
gfortran -shared -fPIC -o scipy/special/amos.so *.o 

gfortran -shared -fPIC -o scipy/special/mach.so *.o

# rm *.o 

mkdir -p $CONDA_PREFIX/lib/scipy/special/
cp scipy/special/amos.so $CONDA_PREFIX/lib/scipy/special/
cp scipy/special/mach.so $CONDA_PREFIX/lib/scipy/special/

# Copy the runtime library

# For Mac
# cp ~/Desktop/lfortran/src/runtime/liblfortran_runtime.dylib $CONDA_PREFIX/lib/

# For Linux
cp ~/Desktop/lfortran/src/runtime/liblfortran_runtime.so $CONDA_PREFIX/lib/

# Build

python dev.py build

# Setup the package

# For Mac
# mkdir -p ./build-install/lib/python3.11/site-packages/scipy/special/
# cp scipy/special/amos.so ./build-install/lib/python3.11/site-packages/scipy/special/


# For linux

mkdir -p ./build-install/lib/python3.11/scipy/special/
cp scipy/special/amos.so ./build-install/lib/python3.11/scipy/special/
cp scipy/special/mach.so ./build-install/lib/python3.11/scipy/special/

# Test

python dev.py test -t scipy.special.tests.test_basic -v

