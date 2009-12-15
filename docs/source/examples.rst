.. ex: set sts=4 ts=4 sw=4 et tw=79:

Examples (Getting Started Guide)
================================

Overview
--------

There are currently two different classes to represent your system and
calculate the entropies. 

DiscreteSystem_
~~~~~~~~~~~~~~~

This is the most general class. It takes two integer arrays of inputs X and Y
and is initiliased in the following way::

    DiscreteSystem(X, X_dims, Y, Y_dims)

X, which usually represents response for neural data, should be of dimensions
``(Xn, T)`` where Xn is the number of variables (cells) and T is the number of
trials. X takes integer values ``0 < X < Xm``. ``X_dims`` is a tuple ``(Xn,
Xm)``. 

Y, which usually represents stimulus for neural data, should be of dimensions
``(Yn, T)`` where Yn is the number of stimulus variables (usually 1) and T is
the number of trials. Y takes integer values ``0 < Y < Ym``. ``Y_dims`` is a
tuple ``(Yn, Ym)``.

.. note::
    In general, X and Y are completely interchangable (due to symmetry of
    mutual information), but it is highly recommended that you choose Y to be
    the space with smaller dimension because of the way the conditional
    distributions are calculated. X and Y can be presented in any order, the
    only requirement is that the trials match up (ie ``X[:,t]`` and ``[Y:,t]``
    both correspond to the values of the system for the same
    trial/observation). 

.. _`DiscreteSystem`: api.html#pyentropy.systems.DiscreteSystem

SortedDiscreteSystem_
~~~~~~~~~~~~~~~~~~~~~

In many cases, the data are available sorted by input (Y value). In this case
the response conditional distributions can be calculated in a much more
efficient way. The class is initiliased in the following way::

    SortedDiscreteSystem(X, X_dims, Ym, Ny)

X is as described above for DiscreteSystem. Ym is the size of of the Y space
finite alphabet (i.e. number of possible Y values) and Ny is an array of length
Ym listing the number of trials available in for each Y value (in the order
that the X values are presented). i.e.  ``X[0:Ny[0]]`` contains trials for Y=0,
``X[Ny[0]:Ny[1]]`` contains trials for Y=1 etc.

.. note::
    All features in pyEntropy are for calculating entropies and mutual 
    information on discrete spaces. This is why the inputs must be integers. 
    These integers are labels for the elements of the discrete space, they 
    have no numerical value. If you're input data is continuous, it must first 
    be quantised using :func:`pyentropy.quantise` to a suitable discrete 
    representation. 

.. _`SortedDiscreteSystem`: api.html#pyentropy.systems.SortedDiscreteSystem

Examples
--------

Here is an example `ipython <http://ipython.scipy.org>`_ session for a simple
example. This example is one of the unit tests. 

.. sourcecode:: ipython

    In [1]: import numpy as np
    In [2]: from pyentropy import DiscreteSystem
    # generate random input
    In [3]: x = np.random.random_integers(0,9,10000)
    # corrupt half of output
    In [4]: y = x.copy()
    In [5]: indx = np.random.permutation(len(x))[:len(x)/2]
    In [6]: y[indx] = np.random.random_integers(0,9,len(x)/2)
    # setup system and calculate entropies
    In [7]: s = DiscreteSystem(x,(1,10), y,(1,10))
    In [8]: s.calculate_entropies(method='plugin', calc=['HX', 'HXY'])
    In [9]: s.I()
    Out[9]: 0.91243680717436071
    In [10]: s.calculate_entropies(method='pt', calc=['HX', 'HXY'])
    In [11]: s.I()
    Out[11]: 0.90659389225876152
    In [12]: s.calculate_entropies(method='qe', calc=['HX', 'HXY'])
    In [13]: s.I()
    Out[13]: 0.90641575663968688
    In [14]: s.calculate_entropies(method='nsb', calc=['HX', 'HXY'])
    In [15]: s.I()
    Out[15]: 0.90307230696261565
    # true value
    In [16]: np.log2(10) - (-9*0.05*np.log2(0.05) - 0.55*np.log2(0.55))
    Out[16]: 0.90268739025051348
    
    
This example simulates a noisy communication channel which randomly corrupts
50% of the observations. 

A more involved example script, which produces the data and plots for Figure 1 
of    
   | Ince, R. A. A., Petersen, R. S., Swan, D. C. and Panzeri, S. (2009)
   | "Python for Information Theoretic Analysis of Neural Data"
   | Frontiers in Neuroinformatics **3**:4 
can be found `here <http://code.google.com/p/pyentropy/wiki/SupplementalData>`_. 
For more information on the entropy values and bias corrections which can be
selected see :ref:`primer`.

