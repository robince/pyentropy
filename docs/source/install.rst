.. ex: set sts=4 ts=4 sw=4 et tw=79:

***************************** 
Installation and Requirements
*****************************


pyEntropy requires a recent version of `Python <http://www.python.org>`_ and
`NumPy <http://www.scipy.org/>`_. Any Python >2.5 and NumPy > 1.3 should be
fine, but the most recent releases are recommend. If you have any problems
please email (or open an issue).

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
    a series of unit tests found in the tests directory.

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



