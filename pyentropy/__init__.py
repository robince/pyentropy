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
A package for calculating bias-corrected estimates of entropy and mutual 
information quantities over discrete spaces.

For more information see the project home page: 
http://code.google.com/p/pyentropy

"""

__author__ = 'Robin Ince'
__version__ = '0.4.0'

from systems import DiscreteSystem, SortedDiscreteSystem
from utils import (prob, decimalise, nsb_entropy, quantise,
                             dec2base, base2dec)

from numpy.testing import Tester
test = Tester().test
