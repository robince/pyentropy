#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Copyright 2008 Robin Ince
import numpy as np
from nose.tools import with_setup, assert_raises
from numpy.testing import *
from pyentropy.utils import *

x = np.arange(3**3)
x2 = np.atleast_2d(x).T
y = np.array([[0, 0, 0, 0],
   [0, 0, 0, 1],
   [0, 0, 0, 2],
   [0, 0, 1, 0],
   [0, 0, 1, 1],
   [0, 0, 1, 2],
   [0, 0, 2, 0],
   [0, 0, 2, 1],
   [0, 0, 2, 2],
   [0, 1, 0, 0],
   [0, 1, 0, 1],
   [0, 1, 0, 2],
   [0, 1, 1, 0],
   [0, 1, 1, 1],
   [0, 1, 1, 2],
   [0, 1, 2, 0],
   [0, 1, 2, 1],
   [0, 1, 2, 2],
   [0, 2, 0, 0],
   [0, 2, 0, 1],
   [0, 2, 0, 2],
   [0, 2, 1, 0],
   [0, 2, 1, 1],
   [0, 2, 1, 2],
   [0, 2, 2, 0],
   [0, 2, 2, 1],
   [0, 2, 2, 2]])

def test_dec2base_1d():
    assert_equal(dec2base(x,3,4),y)

def test_dec2base_2d():
    assert_equal(dec2base(x2,3,4),y)

def test_dec2base_noncol():
    assert_raises(ValueError, dec2base, x2.T, 3, 4)
    
def test_base2dec():
    assert_equal(base2dec(y,3),x)

def test_decimalise():
    assert_equal(decimalise(y.T,4,3),x)

def test_decimalise_error():
    assert_raises(ValueError, decimalise, y, 3, 4)

a1 = np.array([8, 9, 7, 9, 3, 3, 9, 7, 9, 2])
b1 = np.array([0, 0, 1, 2, 0, 0, 0, 2, 1, 4]) / 10.0
a2 = np.array([0, 1, 7, 1, 3, 3, 1, 7, 1, 2])
b2 = np.array([1, 4, 1, 2, 0, 0, 0, 2, 0, 0]) / 10.0

def test_prob_naive():
    assert_equal(prob(a1,10), b1)

def test_prob_naive_missed_responses():
    assert_equal(prob(a2,10), b2)
    
def test_pt_bayescount():
    # values match original bayescount.m file
    for n,r in [(100000, 5.0), (50, 5.0), (30, 6.0),
                (12, 7.0), (10, 8.0), (8, 9.0), (7, 10.0)]:
        yield check_pt_bayes, n, r
        
def check_pt_bayes(n, r):
    assert_equal(pt_bayescount(b1,n),r)
    
if __name__ == '__main__':
    run_module_suite()
