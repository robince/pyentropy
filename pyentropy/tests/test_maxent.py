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
from numpy.testing import *
from pyentropy.maxent import *

# remove cached file
try: 
    os.remove(os.path.join(get_data_dir(),'a_n%im%i.mat'%(3,4)))
except OSError:
    pass
# create from scratch
# need to check both created from scratch and loaded
# to catch any problems with savemat/loadmat round trip
a = AmariSolve(3,4,confirm=False)
# load
a_loaded = AmariSolve(3,4)

# a random distribution
p = np.random.rand(64)
p /= p.sum()

def test_theta_roundtrip():
    assert_array_almost_equal(p, a.p_from_theta(a.theta_from_p(p)))

def test_theta_roundtrip_loaded():
    assert_array_almost_equal(p, a_loaded.p_from_theta(a_loaded.theta_from_p(p)))

# check first order marginals analytic
# this shows numerical solution is accurate
def test_first_order_solve():
    p1a = a.solve(p, 1)
    p1d = order1direct(p, a)
    assert_array_almost_equal(p1a,p1d)

def test_first_order_solve_loaded():
    p1a = a.solve(p, 1)
    p1d = order1direct(p, a_loaded)
    assert_array_almost_equal(p1a,p1d)

if __name__ == '__main__':
    run_module_suite()
