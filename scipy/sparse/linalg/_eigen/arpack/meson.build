# Building ARPACK yields a ton of rank mismatch (scalar and rank-1) warnings
# that cannot be suppressed with a more specific flag.
_suppress_all_warnings = ff.get_supported_arguments('-w')

arpack_module = custom_target('arpack_module',
  output: ['_arpackmodule.c', '_arpack-f2pywrappers.f'],
  input: 'arpack.pyf.src',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

_arpack = py3.extension_module('_arpack',
  arpack_module,
  c_args: numpy_nodepr_api,
  objects: ['arpack.so'],
  link_args: version_link_args,
  dependencies: [lapack, blas, fortranobject_dep, g77_abi_wrappers],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/sparse/linalg/_eigen/arpack'
)

py3.install_sources([
    '__init__.py',
    'arpack.py',
    'ARPACK/COPYING'
  ],
  subdir: 'scipy/sparse/linalg/_eigen/arpack'
)

subdir('tests')
