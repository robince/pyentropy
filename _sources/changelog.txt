.. ex: set sts=4 ts=4 sw=4 et tw=79:

*********
Changelog
*********

0.5.0
-----

* Add Spike Train Analysis Toolkit NSB implementation
* Add Ish-ush estimator
* quantise_discrete for rebinning discrete sequences

0.4.0 - 15/12/09
----------------

* Add HshX and Ish-ush estimator
* Major refactor and cleanup
* ReST docstrings + sphinx documentation
* Unit tests
* Add pyentropy.maxent - finite alphabet maximum entropy solutions

0.3.3 - 7/10/09 
---------------

* Fix serious bug with PiX construction (r54 r57)

0.3.1 - 6/1/09
--------------

* Minor bug for 1D input (r44)
* Add direct Ispike calculation 
* Update docstrings
* Fix nasty bug with shuffled information when used with QE (or if an instance is used more than once with shuffling) (r49)

0.3 - 6/1/09 
------------

* Fix PT correction bug
* Add James-Stein shrinkage estimator
* This version used for plots for Frontiers paper

0.2 - 30/10/08 
--------------

* Make pyentropy a package
* Distutils installation
* Add Pola Decomposition 
* Improve NSB temp file handling
* Add SortedDiscreteSystem class, and refactored common code to base class (BaseSystem)
* Optimise Bayescount function for PT method

0.1 - 21/9/08
-------------

* Initial Release. Basic Python code - `DiscreteSystem` class supports vector inputs, shuffling if necessary for QE.


