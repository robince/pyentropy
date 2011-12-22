.. ex: set sts=4 ts=4 sw=4 et tw=79:

****************
Module Reference
****************

:mod:`pyentropy` -- Core information theoretic functionality
============================================================

.. automodule:: pyentropy.systems

:class:`BaseSystem`
-------------------------------------

.. autoclass:: pyentropy.systems.BaseSystem
   :members:

:class:`DiscreteSystem`
--------------------------------

.. autoclass:: pyentropy.systems.DiscreteSystem
   :show-inheritance:
   :members: __init__

:class:`SortedDiscreteSystem`
---------------------------------------

.. autoclass:: pyentropy.systems.SortedDiscreteSystem
   :show-inheritance:
   :members: __init__

Utility Functions
-----------------

.. autofunction:: pyentropy.base2dec

.. autofunction:: pyentropy.dec2base

.. autofunction:: pyentropy.decimalise

.. autofunction:: pyentropy.nsb_entropy

.. autofunction:: pyentropy.prob

.. autofunction:: pyentropy.quantise

.. autofunction:: pyentropy.quantise_discrete



:mod:`pyentropy.maxent` -- Finite Alphabet Maximum-Entropy Solutions
====================================================================

.. automodule:: pyentropy.maxent

.. autoclass:: pyentropy.maxent.AmariSolve
   :members: __init__, solve, theta_from_p, eta_from_p, p_from_theta

.. autofunction:: pyentropy.maxent.get_config_file

.. autofunction:: pyentropy.maxent.get_data_dir

:mod:`pyentropy.statk` -- Python wrapper of `STA Toolkit`_ functions
====================================================================

.. automodule:: pyentropy.statk

.. autofunction:: pyentropy.statk.nsb_entropy

.. autofunction:: pyentropy.statk.bub_entropy

.. _STA Toolkit: http://neuroanalysis.org/neuroanalysis/goto.do?page=.repository.toolkit_home
