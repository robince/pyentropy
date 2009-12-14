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
"""
Entropy and Information Estimates
=================================

 This module provides techniques for estimates of information theoretic 
 quantities.

Classes
-------

  DiscreteSystem :
      Class to sample and hold probabilities of a general discrete system.
  SortedDiscreteSystem:
      Class to sample and hold probabilites for a system where the input
      output mapping is available already sorted by output.

"""
__author__ = 'Robin Ince'
__version__ = '0.4.0dev'

from systems import DiscreteSystem, SortedDiscreteSystem
from utils import (prob, decimalise, nsb_entropy, quantise,
                             dec2base, base2dec)

from numpy.testing import Tester
test = Tester().test
