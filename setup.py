import sys, os
import pyentropy

# BEFORE importing disutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

from distutils.core import setup
from distutils.extension import Extension
from distutils.command.build_ext import build_ext
from distutils.errors import CCompilerError, DistutilsExecError

extension_build_failed = False

def ext_failed_warning(name):
    print ('*'*70+'\n')*3

    print """WARNING: The %s extension module could not be 
compiled. pyEntropy should run, but the features
present in that file will not be available.

Above is the ouput showing how the compilation
failed."""%name

    if sys.platform == 'win32':
        print

        print """I see you are using Windows. The default
compiler for this platform is the Microsoft Visual
Studio C compiler. However, a free alternative
compiler called mingw can be used instead."""

    print
    print ('*'*70+'\n')*3
    global extension_build_failed
    extension_build_failed = True

try:
    from gsl_dist.gsl_Extension import gsl_Extension
except DistutilsExecError:
    ext_failed_warning('gsl-based')

exts = []
wrap_sources = ['hist_c.c', 'sort_c.c', 'gen_c.c', 'entropy_c.c',
                      'entropy_nsb_c.cpp', 'entropy_bub_c.c', 'wrap.c']
statk_wrap_sources = [os.path.join('pyentropy','statk',x) for x in wrap_sources]
try: 
    statk_wrap = gsl_Extension("statk.wrap",
                           sources = statk_wrap_sources,
                           gsl_min_version=(1,),
                           python_min_version=(2,5)
                           )
    exts.append(statk_wrap)                           
except:
    pass

class build_ext_allow_fail( build_ext ):
    # This class allows C extension building to fail.
    # Taken from visionegg (LGPL)
    # http://github.com/visionegg/visionegg/blob/master/setup.py
    # http://www.visionegg.org/
    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except CCompilerError, x:
            ext_failed_warning(ext.name)


setup(name='pyentropy',
      version=pyentropy.__version__,
      description='Entropy and Information Theoretic Estimates',
      author=pyentropy.__author__,
      author_email='pyentropy@robince.net',
      url='http://code.google.com/p/pyentropy',
      packages=['pyentropy','pyentropy.tests','pyentropy.statk'],
      ext_package='pyentropy',
      ext_modules=exts,
      cmdclass={'build_ext':build_ext_allow_fail}
      )

if extension_build_failed:
    print ('*'*70+'\n')*3

    print """WARNING: Building of some extensions failed. Please
see the messages above for details.\n"""

    print ('*'*70+'\n')*3

