#
# python3 setup_lifo.py install --home=.
#
from distutils.core import setup, Extension
from distutils.sysconfig import get_python_lib

def main():
  import os
  incdir = os.path.join(get_python_lib(plat_specific=1), "numpy/core/include")
  mach = os.path.basename(os.getcwd())
  setup(name="Lifo",
          version="1.0.0",
          description="Python interface for the lifo C library function",
          author="Alain Hebert",
          author_email="alain.hebert@polymtl.ca",
          ext_modules=[Extension("lifo", ["lifomodule.c"],
                     include_dirs=["../../../Ganlib/src",incdir],
                     library_dirs=["../../lib/"+mach],
                     libraries=["Ganlib"] ) ])

if __name__ == "__main__":
    main()
