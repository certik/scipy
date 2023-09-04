# TODO: Add flags for 64 bit ints

interpnd = py3.extension_module('interpnd',
  spt_cython_gen.process('interpnd.pyx'),
  c_args: [Wno_discarded_qualifiers, cython_c_args],
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/interpolate'
)

py3.extension_module('_ppoly',
  linalg_cython_gen.process('_ppoly.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/interpolate'
)

py3.extension_module('_rgi_cython',
  linalg_cython_gen.process('_rgi_cython.pyx'),
  c_args: cython_c_args,
  dependencies: np_dep,
  install: true,
  subdir: 'scipy/interpolate'
)

_bspl = py3.extension_module('_bspl',
  cython_gen.process('_bspl.pyx'),
  c_args: cython_c_args,
  include_directories: 'src/',
  dependencies: np_dep,
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/interpolate'
)

# TODO: Add flags for 64 bit ints
_fitpack = py3.extension_module('_fitpack',
  ['src/_fitpackmodule.c'],
  c_args: numpy_nodepr_api,
  objects: 'fitpack.so',
  include_directories: 'src/',
  dependencies: np_dep,
  link_args: version_link_args,
  link_language: 'fortran',
  install: true,
  subdir: 'scipy/interpolate'
)

dfitpack_module = custom_target('dfitpack_module',
  output: ['dfitpack-f2pywrappers.f', 'dfitpackmodule.c'],
  input: 'src/fitpack.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

# TODO: Add flags for 64 bit ints
dfitpack = py3.extension_module('dfitpack',
  dfitpack_module,
  c_args: [numpy_nodepr_api, Wno_unused_variable],
  link_args: version_link_args,
  dependencies: [lapack, fortranobject_dep],
  objects: 'fitpack.so',
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/interpolate'
)

if use_pythran
  _rbfinterp_pythran = custom_target('_rbfinterp_pythran',
    output: ['_rbfinterp_pythran.cpp'],
    input: '_rbfinterp_pythran.py',
    command: [pythran, '-E', '@INPUT@', '-o', '@OUTDIR@/_rbfinterp_pythran.cpp']
  )

  _rbfinterp_pythran = py3.extension_module('_rbfinterp_pythran',
    _rbfinterp_pythran,
    cpp_args: cpp_args_pythran,
    dependencies: [pythran_dep, np_dep],
    link_args: version_link_args,
    install: true,
    subdir: 'scipy/interpolate'
  )
else
  py3.install_sources(
    ['_rbfinterp_pythran.py'],
    subdir: 'scipy/interpolate'
  )
endif

py3.install_sources([
    '__init__.py',
    '_bsplines.py',
    '_cubic.py',
    '_fitpack_impl.py',
    '_fitpack_py.py',
    '_fitpack2.py',
    '_interpnd_info.py',
    '_interpolate.py',
    '_ndgriddata.py',
    '_pade.py',
    '_polyint.py',
    '_rbf.py',
    '_rbfinterp.py',
    '_rgi.py',
    '_ndbspline.py',
    'fitpack.py',
    'fitpack2.py',
    'interpolate.py',
    'ndgriddata.py',
    'polyint.py',
    'rbf.py'
  ],
  subdir: 'scipy/interpolate'
)

subdir('tests')
