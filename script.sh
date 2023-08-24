lfortran -c scipy/special/specfun/specfun.f -o scipy/special/specfun.o --fixed-form --implicit-typing --implicit-interface --generate-object-code  --use-loop-variable-after-loop --rtlib
clang -shared -fPIC -o scipy/special/specfun.so scipy/special/specfun.o -L/Users/pranavchiku/repos/lfortran/src/runtime -llfortran_runtime
cp scipy/special/specfun.so $CONDA_PREFIX/lib/scipy/special/
cp scipy/special/specfun.so ./build-install/lib/python3.10/site-packages/scipy/special/
python dev.py test -t scipy.special.tests.test_hyp2f1 -v
