.. ex: set sts=4 ts=4 sw=4 et tw=79:

.. _primer: 

***********************************************
Entropy, Information and Bias Correction Primer
***********************************************

Entropy and Mutual Information
==============================

Information Breakdown
---------------------

Bias
====

Origins of sampling bias
------------------------

Bias correction methods
-----------------------

pyEntropy
=========

Bias Corrections
----------------

pyEntropy currently implements the following bias correction methods, specified
with `method` argument to
:func:`pyentropy.systems.BaseSystem.calculate_entropies`:

``plugin``
    The standard *plugin* or *naive* entropy estimator.
``pt``
    Panzeri-Treves bias correction [PT96]_. This is the Miller-Madow analytic
    correction, where the number of non-zero responses is estimated
    using a Bayesian procedure rather than the naive count, yielding much
    better results. 
``qe``
    Quadratic extrapolation [Strong98]_. See above.
``nsb``
    Nemenman-Shafee-Bialek method [NSB02]_. This uses a C++ version of the
    algorithm from the `Spike Train Analysis Toolkit 
    <http://neuroanalysis.org/toolkit/index.html>`_. 
``nsb-ext``
    Nemenman-Shafee-Bialek method [NSB02]_. This uses the external
    ``nsb-entropy`` package from http://nsb-entropy.sourceforge.net/ . This
    should be installed on your path.

Which method should I use?
~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a difficult question and is to some extent subjective, being a trade
off between bias, variance and computational requirements. It is complicated by
the fact that these factors are also affected by the statistics of the system
being studied. Generally is you have a single variable X space the best results
are obtained with the NSB estimator, but it is much slower to compute than the
other methods. If you have a multi-variable X space, the best results are
achieved with the shuffled estimator, corrected with either PT or QE. The best
thing to do, if possible, is to test the methods with simulated data with
similar statistics. Future versions of pyEntropy will include helper functions
to make this easier.

Entropy Values
--------------

pyEntropy currently implements the following entropy values, specified in
`calc` argument to :func:`pyentropy.systems.BaseSystem.calculate_entropies`:

``HX``
    :math:`H(\mathbf{X})` -- unconditional entropy of X.
``HY``
    :math:`H(Y)` -- unconditional entropy of Y.
``HshX``
    :math:`H_{sh}(\mathbf{X})` -- shuffle-decorrelated independent unconditional    entropy of X. (X components are shuffled). 
``SiHXi``
    :math:`\sum_{i=1}^{Xm} H(X_{i})` -- direct-decorrelated independent 
    unconditional entropy X. Summing the individual entropies corresponds 
    to constructing the analytic independent joint distribution (product 
    of individual distributions).
``HXY``
    :math:`H(\mathbf{X}|Y)` -- entropy of X conditional on Y.
``HiXY``
    :math:`H_{ind}(\mathbf{X}|Y) = \sum_{i=1}^{Xm} H(X_{i}|Y)` --
    direct-decorrelated independent conditional entropy. 
``HshXY``
    :math:`H_{sh}(\mathbf{X}|Y)` -- shuffle-decorrelated independent
    conditional entropy. X components are shuffled for each response Y (but not
    between responses). 
``HiX``
    :math:`H_{ind}(\mathbf{X})` -- unconditional direct-conditionally-decorrelated entropy of X. Entropy of :math:`P_{ind}(X) = \sum_{y \in Y}
    P(y) \prod_{i=1}^{Xm} P(X_{i}|y)`
``ChiXY``
    :math:`\chi (\mathbf{X})` -- cross entropy between :math:`P_{ind}(X)` and 
    :math:`P(X)`
