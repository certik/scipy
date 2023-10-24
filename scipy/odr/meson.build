py3.extension_module('__odrpack',
  '__odrpack.c',
  objects: 'odrpack.so',
  link_args: version_link_args,
  c_args: numpy_nodepr_api,
  dependencies: [blas, np_dep],
  install: true,
  link_language: 'fortran',
  subdir: 'scipy/odr'
)

python_sources = [
  '__init__.py',
  '_add_newdocs.py',
  '_models.py',
  '_odrpack.py',
  'models.py',
  'odrpack.py'
]

py3.install_sources(
  python_sources,
  subdir: 'scipy/odr'
)

subdir('tests')
