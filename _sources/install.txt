.. ex: set sts=4 ts=4 sw=4 et tw=79:

***************************** 
Installation and Requirements
*****************************


pyEntropy requires a recent version of `Python <http://www.python.org>`_ and
`NumPy <http://www.scipy.org/>`_. Any Python >2.5 (but < 3.0) and NumPy > 1.3
should be fine, but the most recent releases are recommend. If you have any
problems please email (or open an issue). From version 0.5.0, a C++ NSB
implementation (from the `Spike Train Analysis Toolkit
<http://neuroanalysis.org/toolkit/index.html>`_) is included which requires a
compiler, and also the `Gnu Scientific Library (GSL)
<http://www.gnu.org/software/gsl/>`_. On Linux, you should be able to install
these through your package manager. E.g. for Ubuntu::

    sudo aptitude install build-essential python-dev libgsl0-dev

should provide all the requirements. On a Mac, XCode and (with MacPorts
installed)::

    sudo port install gsl

should be enough. On Posix platforms (Mac/Linux) you should have the
``gsl-config`` program in your path. For Windows, I am currently investigating
the best way, and will hopefully provide binaries.

.. note:: 

    The :mod:`pyentropy.statk.wrap` module is optional. I intend to keep
    pyEntropy installable and usable as a pure Python package with minimal
    dependencies. If the build fails for any reason, warnings will be
    displayed, but the pyEntropy package should still be installable. The
    ``nsb`` entropy method will be unavailable but if you have the
    ``nsb-entropy`` program from http://nsb-entropy.sourceforge.net/ on your
    path you should still be able to use the ``nsb-ext`` method.

The :mod:`pyentropy.maxent` module also requires 
`SciPy <http://www.scipy.org>`_.

For windows, running the installer should be all that is needed. On other
platforms, uncompress the archive and run the following command::

    python setup.py install
    
This will install the package to the appropriate location of your Python
installation. This package uses distutils; more details are available in the
`python documentation <http://docs.python.org/install/index.html>`_. Note that
from Python 2.6 you can run::
    
    python setup.py install --user 

to install to ``.local`` directory in the home directory. This is automatically
added to PYTHONPATH and is the recommended way to install for a single user on
UNIX systems (without root access).

.. note::
    You can test your installation by running ``pyentropy.test()``. This runs
    a series of unit tests found in the tests directory. Since the tests of the
    NSB method can take a long time (around 2 minutes on a fast computer) they
    are not run by default, but can be included by running
    ``pyentropy.test('full')``. 

:mod:`pyentropy.maxent` stores generated transformation matrices for a given
set of parameters are to disk. These files can get large so you should be aware
that they are there. The default location for the cache is a ``.pyentropy``
(``_pyentropy`` on windows) directory in the users home directory. To override
this and use a custom location (for example to share the folder between users)
you can put a configuration file ``.pyentropy.cfg`` (``pyentropy.cfg`` on
windows) file in the home directory with the following format::

    [maxent]
    cache_dir = /path/to/cache
    
:func:`pyentropy.maxent.get_config_file()` will show where it is looking for the config
file.



