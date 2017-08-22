#    This file is part of pyEntropy
#
#    pyEntropy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    pyEntropy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with pyEntropy. If not, see <http://www.gnu.org/licenses/>.
#
#    Copyright 2009, 2010 Robin Ince

from copy import copy
import numpy as np
from numpy.testing import *
from nose.tools import with_setup
from pyentropy import DiscreteSystem, SortedDiscreteSystem
import subprocess, os

# TODO: test ChiXY1 HXY1 for binary data (Adelman Ispike)
# TODO: test running more than once on an instance (to catch eg shuffling bug)
# TODO: explicitly test shuffling?
# TODO: test performance of bias corrections?
# TODO: test decomposition helpers

def setup():
    global allcalc
    # all entropy values
    allcalc = ['HX','HY','HXY','SiHXi','HiX','HshX','HiXY',
               'HshXY','ChiX']

def statk_available():
    try:
        import pyentropy.statk.wrap
    except ImportError:
        return False
    else:
        return True

def nsb_ext_available():
    try:
        subprocess.call(["nsb-entropy"],stdout=open(os.devnull))
    except OSError:
        return False
    else:
        return True
        
#
# test DiscreteSystem with simple 1D input, output
#
def setup_1d():
    global x, y, HX, HXY, alltrue
    # simple channel which corrupts 50% randomly
    x = np.random.randint(0,9+1,100000)
    y = x.copy()
    indx = np.random.permutation(len(x))[:len(x)/2]
    y[indx] = np.random.randint(0,9+1,len(x)/2)
    # analytic results
    HX = np.log2(10)
    HXY = -9*0.05*np.log2(0.05) - 0.55*np.log2(0.55)
    # same order as allcalc
    alltrue = np.array([HX, HX, HXY, HX, HX, HX, HXY, HXY, HX])

def teardown_1d():
    global x, y, HX, HXY, alltrue
    del x, y, HX, HXY, alltrue

@with_setup(setup_1d, teardown_1d)
def do_1d_check(method, qe_method):
    xc = x.copy()
    yc = y.copy()
    s = DiscreteSystem(xc,(1,10),yc,(1,10))
    # calculate all entropies
    s.calculate_entropies(method=method, calc=allcalc, qe_method=qe_method)
    # check output assinged
    assert_(s.H == getattr(s,'H_%s'%method.replace('-','')))
    v = np.array([s.H[k] for k in allcalc])
    assert_array_almost_equal(v, alltrue, decimal=2)
    # check didn't do something nasty to inputs
    assert_array_equal(x, xc)
    assert_array_equal(y, yc)

def test_1d_plugin():
    yield do_1d_check, 'plugin', None
    
def test_1d_pt():
    yield do_1d_check, 'pt', None

@dec.skipif(not statk_available(), "STATK wrapper not available")
@dec.slow
def test_1d_nsb():
    yield do_1d_check, 'nsb', None

@dec.skipif(not nsb_ext_available(), "nsb-entropy binary not found on path")
@dec.slow
def test_1d_nsbext():
    yield do_1d_check, 'nsb-ext', None

@dec.skipif(not statk_available(), "STATK wrapper not available")
@dec.slow
def test_1d_bub():
    yield do_1d_check, 'bub', None

def test_1d_qe():
    yield do_1d_check, 'qe', 'plugin'

def test_1d_qe_pt():
    yield do_1d_check, 'qe', 'pt'



#
# test SortedDiscreteSystem with simple 1D input, output
#
def setup_1d_sorted():
    global x, y, Ny
    setup_1d()
    # convert to sorted system format
    xs = np.zeros_like(x)
    Ny = np.zeros(10)
    start = 0
    for i in range(10):
        oce = x[y==i]
        Ny[i] = len(oce)
        end = start + Ny[i]
        xs[int(start):int(end)] = oce
        start = end
    x = xs

def teardown_1d_sorted():
    global x, y, Ny, HX, HXY, alltrue
    del x, y, Ny, HX, HXY, alltrue

@with_setup(setup_1d_sorted, teardown_1d_sorted)
def do_1d_check_sorted(method, qe_method):
    xc = x.copy()
    yc = y.copy()
    s = SortedDiscreteSystem(xc,(1,10),10, Ny)
    # calculate all entropies
    s.calculate_entropies(method=method, calc=allcalc, qe_method=qe_method)
    # check output assinged
    assert_(s.H == getattr(s,'H_%s'%method.replace('-','')))
    v = np.array([s.H[k] for k in allcalc])
    assert_array_almost_equal(v, alltrue, decimal=2)
    # check didn't do something nasty to inputs
    assert_array_equal(x, xc)
    assert_array_equal(y, yc)

def test_1d_plugin_sorted():
    yield do_1d_check_sorted, 'plugin', None
    
def test_1d_pt_sorted():
    yield do_1d_check_sorted, 'pt', None

@dec.skipif(not statk_available(), "STATK NSB wrapper not available")
@dec.slow
def test_1d_nsb_sorted():
    yield do_1d_check_sorted, 'nsb', None

@dec.skipif(not nsb_ext_available(), "nsb-entropy binary not found on path")
@dec.slow
def test_1d_nsbext_sorted():
    yield do_1d_check_sorted, 'nsb-ext', None

@dec.skipif(not statk_available(), "STATK wrapper not available")
@dec.slow
def test_1d_bub_sorted():
    yield do_1d_check_sorted, 'bub', None

def test_1d_qe_sorted():
    yield do_1d_check_sorted, 'qe', 'plugin'

def test_1d_qe_pt_sorted():
    yield do_1d_check_sorted, 'qe', 'pt'

#
# toy system to check decomposition, PiX construction etc.
# 

def setup_toy1():
    global x, y, Ny, toycalc, toytrue
    x = np.array([[0,1,1], [1,1,2], [1,1,1], [0,1,1], [0,2,1], [0,0,0]]).T
    y = np.array([0,0,0,1,1,1])
    Ny = np.array([3,3])
    toycalc = ['HX','HXY','HiXY','ChiX', 'HiX', 'SiHXi']
    # true values checked against ibtb
    toytrue = np.array([2.2516291673878226,
                        1.5849625007211561,
                        2.1699250014423122,
                        2.8365916681089796,
                        2.9477027792200903,
                        3.4215541688301352])

def teardown_toy1():
    global x, y, Ny, toycalc
    del x, y, Ny, toycalc

@with_setup(setup_toy1, teardown_toy1)
def test_toy1():
    s = DiscreteSystem(x,(3,3),y,(1,2))
    s.calculate_entropies(method='plugin', calc=toycalc)
    v = np.array([s.H[t] for t in toycalc])
    assert_array_almost_equal(v, toytrue)

@with_setup(setup_toy1, teardown_toy1)
def test_toy1_sorted():
    s = SortedDiscreteSystem(x,(3,3),2,Ny)
    s.calculate_entropies(method='plugin', calc=toycalc)
    v = np.array([s.H[t] for t in toycalc])
    assert_array_almost_equal(v, toytrue)

    
if __name__ == '__main__':
    run_module_suite()
