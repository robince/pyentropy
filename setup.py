from distutils.core import setup
import pyentropy

setup(name='pyentropy',
      version=pyentropy.__version__,
      description='Entropy and Information Theoretic Estimates',
      author=pyentropy.__author__,
      author_email='pyentropy@robince.net',
      url='http://code.google.com/p/pyentropy',
      packages=['pyentropy','pyentropy.tests']
      )
