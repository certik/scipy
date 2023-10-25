lsoda_src = [
  'odepack/blkdta000.f',
  'odepack/bnorm.f',
  'odepack/cfode.f',
  'odepack/ewset.f',
  'odepack/fnorm.f',
  'odepack/intdy.f',
  'odepack/lsoda.f',
  'odepack/prja.f',
  'odepack/solsy.f',
  'odepack/srcma.f',
  'odepack/stoda.f',
  'odepack/vmnorm.f',
  'odepack/xerrwv.f',
  'odepack/xsetf.f',
  'odepack/xsetun.f'
]

vode_src = [
  'odepack/vode.f',
  'odepack/zvode.f'
]

dop_src = [
  'dop/dop853.f',
  'dop/dopri5.f'
]

quadpack_test_src = [
  'tests/_test_multivariate.c'
]

lsoda_lib = static_library('lsoda_lib',
  lsoda_src,
  fortran_args: fortran_ignore_warnings
)

vode_lib = static_library('vode_lib',
  vode_src,
  fortran_args: fortran_ignore_warnings
)

dop_lib = static_library('dop_lib',
  dop_src,
  fortran_args: fortran_ignore_warnings
)

py3.extension_module('_quadpack',
  ['__quadpack.h', '_quadpackmodule.c'],
  c_args: numpy_nodepr_api,
  fortran_args: fortran_ignore_warnings,
  include_directories: ['../_lib/src'],
  objects: ['quadpack.so', 'integrate_mach.so'],
  link_args: version_link_args,
  dependencies: [lapack, np_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

py3.extension_module('_odepack',
  '_odepackmodule.c',
  c_args: numpy_nodepr_api,
  fortran_args: fortran_ignore_warnings,
  objects: 'integrate_mach.so',
  link_with: [lsoda_lib],
  link_args: version_link_args,
  dependencies: [lapack, np_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

vode_module = custom_target('vode_module',
  output: ['_vode-f2pywrappers.f', '_vodemodule.c'],
  input: 'vode.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('_vode',
  vode_module,
  link_with: [vode_lib],
  c_args: [numpy_nodepr_api, Wno_unused_variable],
  link_args: version_link_args,
  dependencies: [lapack, blas, fortranobject_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

lsoda_module = custom_target('lsoda_module',
  output: ['_lsoda-f2pywrappers.f', '_lsodamodule.c'],
  input: 'lsoda.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('_lsoda',
  lsoda_module,
  objects: 'integrate_mach.so',
  link_with: [lsoda_lib],
  c_args: [numpy_nodepr_api, Wno_unused_variable],
  dependencies: [lapack, fortranobject_dep],
  link_args: version_link_args,
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

_dop_module = custom_target('_dop_module',
  output: ['_dop-f2pywrappers.f', '_dopmodule.c'],
  input: 'dop.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('_dop',
  _dop_module,
  link_with: [dop_lib],
  c_args: [numpy_nodepr_api, Wno_unused_variable],
  dependencies: [fortranobject_dep],
  link_args: version_link_args,
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

py3.extension_module('_test_multivariate',
  [quadpack_test_src],
  link_args: version_link_args,
  install: true,
  subdir: 'scipy/integrate'
)

_test_odeint_banded_module = custom_target('_test_odeint_banded_module',
  output: ['_test_odeint_bandedmodule.c', '_test_odeint_banded-f2pywrappers.f'],
  input: 'tests/banded5x5.pyf',
  command: [py3, generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@']
)

py3.extension_module('_test_odeint_banded',
  ['tests/banded5x5.f', _test_odeint_banded_module],
  c_args: numpy_nodepr_api,
  objects: 'integrate_mach.so',
  link_with: [lsoda_lib],
  fortran_args: _fflag_Wno_unused_dummy_argument,
  link_args: version_link_args,
  dependencies: [lapack, fortranobject_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/integrate'
)

subdir('_ivp')
subdir('tests')

py3.install_sources([
    '__init__.py',
    '_bvp.py',
    '_ode.py',
    '_quad_vec.py',
    '_quadrature.py',
    '_odepack_py.py',
    '_quadpack_py.py',
    '_tanhsinh.py',
    'dop.py',
    'lsoda.py',
    'odepack.py',
    'quadpack.py',
    'vode.py',
  ],
  subdir: 'scipy/integrate'
)
