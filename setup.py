# setup.py

from distutils.core import setup, Extension
import distutils.sysconfig
import numpy

cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:        
        cfg_vars[key] = value.replace("-O2", "")

for key, value in cfg_vars.items():
    if type(value) == str:        
        cfg_vars[key] = value.replace("-DNDEBUG", "")

cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:        
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")


        
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

rt_module = Extension('_rt', libraries=['CGAL','gmp','gomp'],
								extra_compile_args = ["-g","-ffast-math","-fopenmp","-std=c++11","-fopenmp"],
							    sources=['integrate.cxx','integrateQ1.cxx','integrateQ0.cxx','edges.cxx','laguerre.cxx','rt2.cxx','rt.i'],
								include_dirs = [numpy_include],swig_opts=['-c++','-py3'])

setup(name='rt', ext_modules=[rt_module], py_modules=["rt"])
