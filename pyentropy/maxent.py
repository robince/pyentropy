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
Module for computing finite-alphabet maximum entropy solutions using a 
coordinate transform method

For details of the method see:

    Ince, R. A. A., Petersen, R. S., Swan, D. C., Panzeri, S., 2009
    "Python for Information Theoretic Analysis of Neural Data", 
    Frontiers in Neuroinformatics 3:4 doi:10.3389/neuro.11.004.2009
    http://www.frontiersin.org/neuroinformatics/paper/10.3389/neuro.11/004.2009/
    
If you use this code in a published work, please cite the above paper.

The generated transformation matrices for a given set of parameters are 
stored to disk. The default location for the cache is a ``.pyentropy`` 
(``_pyentropy`` on windows) directory in the users home directory. To 
override this and use a custom location (for example to share the folder 
between users) you can put a configuration file ``.pyentropy.cfg`` 
(``pyentropy.cfg`` on windows) file in the home directory with the 
following format::

    [maxent]
    cache_dir = /path/to/cache
    
:func:`pyentropy.maxent.get_config_file()` will show where it is looking for the config
file.

The probability vectors for a finite-alphabet space of ``n`` variables with
``m`` possible values is a length ``m**n-1`` vector ordered such that the 
value of the index is equal to the decimal value of the input state 
represented, when interpreted as a base m, length n word. eg for n=3,m=3::

    P[0] = P(0,0,0)
    P[1] = P(0,0,1)
    P[2] = P(0,0,2)
    P[3] = P(0,1,0)
    P[4] = P(0,1,1) etc.

This allows efficient vectorised conversion between probability index and 
response word using base2dec, dec2base. The output is in the same format.

"""
import time
import os
import sys
import cPickle
import numpy as np
import scipy as sp
import scipy.io as sio
import scipy.sparse as sparse
import scipy.optimize as opt
# umfpack disabled due to bug in scipy
# http://mail.scipy.org/pipermail/scipy-user/2009-December/023625.html
#try:
    #import scikits.umfpack as um
    #HAS_UMFPACK = True
#except:
    #HAS_UMFPACK = False
HAS_UMFPACK = False
from scipy.sparse.linalg import spsolve, use_solver
use_solver(useUmfpack=False)
from utils import dec2base, base2dec
import ConfigParser

def get_config_file():
    """Get the location and name of the config file for specifying
    the data cache dir. You can call this to find out where to put your
    config.

    """
    if sys.platform.startswith('win'):
        cfname = '~/pyentropy.cfg'
    else:
        cfname = '~/.pyentropy.cfg'
    return os.path.expanduser(cfname)

def get_data_dir():
    """Get the data cache dir to use to load and save precomputed matrices"""
    # default values
    if sys.platform.startswith('win'):
        dirname = '~/_pyentropy'
    else:
        dirname = '~/.pyentropy'
    # try to load user override
    config = ConfigParser.RawConfigParser()
    cf = config.read(get_config_file())
    try:
        data_dir = os.path.expanduser(config.get('maxent','cache_dir'))
    except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
        data_dir = os.path.expanduser(dirname)

    # check directory exists
    if not os.path.isdir(data_dir):
        try:
            os.mkdir(data_dir)
        except:
            print "ERROR: could not create data dir. Please check your " + \
                  "configuration."
            raise
    return data_dir

#
# AmariSolve class
#
class AmariSolve:
    """A class for computing maximum-entropy solutions.

    When the class is initiliased the coordinate transform matrices are loaded
    from disk, if available, or generated.

    See module docstring for location of cache directory.
    
    An instance then exposes a solve method which returns the maximum entropy
    distribution preserving marginal constraints of the input probability 
    vector up to a given order k. 

    This class computed the full transformation matrix and so can compute 
    solutions for any order.
   
    """

    def __init__(self, n, m, filename='a_', local=False, confirm=True):
        """Setup transformation matrix for given parameter set.

        If existing matrix file is found, load the (sparse) transformation
        matrix A, otherwise generate it.

        :Parameters:
          n : int
            number of variables in the system
          m : int
            size of finite alphabet (number of symbols)
          filename : {str, None}, optional
            filename to load/save (designed to be used by derived classes).
          local : {False, True}, optional 
            If True, then store/load arrays from 'data/' directory in 
            current working directory. Otherwise use the package data dir
            (default ~/.pyentropy or ~/_pyentropy (windows))
            Can be overridden through ~/.pyentropy.cfg or ~/pyentropy.cfg 
            (windows)
          confirm : {True, False}, optional
            Whether to prompt for confirmation before generating matrix

        """

        #if np.mod(m,2) != 1:
         #   raise ValueError, "m must be odd"

        try: 
            k = self.k
        except AttributeError:
            self.k = n
            
        self.n = n
        self.m = m
        self.l = (m-1)/2
        # full dimension of probability space
        self.fdim = m**n
        # dimension of arrays (-1 dof)
        self.dim = self.fdim - 1

        filename = filename + "n%im%i"%(n,m) 
        if local:
            self.filename = os.path.join(os.getcwd(), 'data', filename)
        else:
            self.filename = os.path.join(get_data_dir(), filename)

        # if file exists load (matrix A)
        # must be running in correct directory
        if os.path.exists(self.filename+'.mat'):
            loaddict = sio.loadmat(self.filename+'.mat')
            self.A = loaddict['A'].tocsc()
            self.order_idx = loaddict['order_idx'].squeeze()
        elif confirm:
            inkey = raw_input("Existing .mat file not found..." +
                              "Generate matrix? (y/n)")
            if inkey == 'y':
                # else call matrix generation function (and save)
                self._generate_matrix()
            else:
                print "File not found and generation aborted..."
                print "Do not use this class instance."
                return None
        else:
            # just generate it without confirmation
            self._generate_matrix()

        self.B = self.A.T
        # umfpack factorisation of matrix
        if HAS_UMFPACK:
            self._umfpack()

        return None

    def _umfpack(self):
        self.umf = um.UmfpackContext()
        self.umf.numeric(self.B)

    def _calculate_orders(self):
        k = self.k
        n = self.n
        m = self.m
        dim = self.dim
        
        # Calculate the length of each order
        self.order_idx       = np.zeros(n+2, dtype=int) 
        self.order_length    = np.zeros(n+1, dtype=int)
        self.row_counter     = 0

        for ordi in xrange(n+1):    
            self.order_length[ordi] = (sp.misc.comb(n, ordi+1, exact=1) *
                                        ((m-1)**(ordi+1)))
            self.order_idx[ordi] = self.row_counter
            self.row_counter += self.order_length[ordi]

        self.order_idx[n+1] = dim+1

        # Calculate nnz for A
        # not needed for lil sparse format
        x = (m*np.ones(n))**np.arange(n-1,-1,-1)
        x = x[:k]
        y = self.order_length[:k]
        self.Annz = np.sum(x*y.T)
        
    def _generate_matrix(self):
        """Generate A matrix if required"""
        k = self.k
        n = self.n
        m = self.m
        dim = self.dim

        self._calculate_orders()

        self.A = sparse.dok_matrix((self.order_idx[k],dim))

        self.row_counter = 0
        for ordi in xrange(k):
            self.nterms = m**(n - (ordi+1))
            self.terms = dec2base(np.c_[0:self.nterms,], m, n-(ordi+1))
            self._recloop((ordi+1), 1, [], [], n, m)
            print "Order " + str(ordi+1) + " complete. Time: " + time.ctime()

        # save matrix to file
        self.A = self.A.tocsc()
        savedict = {'A':self.A, 'order_idx':self.order_idx}
        sio.savemat(self.filename, savedict)

    def _recloop(self, order, depth, alpha, pos, n, m, blocksize=None):
        terms = self.terms
        A = self.A
        if not blocksize:
            blocksize = self.nterms

        # starting point for position loop
        if len(pos)==0:
            pos_start = 0
        else:
            pos_start = pos[-1] + 1

        # loop over alphabet
        for ai in xrange(1, m):
            alpha_new = list(alpha)
            alpha_new.append(ai)

            # loop over position
            for pi in xrange(pos_start, (n-(order-depth))):
                pos_new = list(pos)
                pos_new.append(pi)

                # add columns?
                if depth == order:
                    # special case for highest order
                    # (can't insert columns into empty terms array)
                    if order==n:
                        cols = base2dec(np.atleast_2d(alpha_new),m)[0]-1
                        A[self.row_counter, cols] = 1
                    else:    
                        # add columns (insert and add to sparse)
                        ins = np.tile(alpha_new,(blocksize,1))
                        temp = terms
                        for coli in xrange(order):
                            temp = inscol(temp, np.array(ins[:,coli],ndmin=2).T, pos_new[coli])

                        cols = (base2dec(temp,m)-1).tolist()

                        A[self.row_counter, cols] = 1;

                    self.row_counter += 1
                else:
                    self._recloop(order, depth+1, alpha_new, pos_new, n, m, blocksize=blocksize)

    def solve(self,Pr,k,eta_given=False,ic_offset=-0.01, **kwargs):
        """Find maxent distribution for a given order k
        
        :Parameters:
          Pr : (fdim,)
           probability distribution vector
          k : int
            Order of interest (marginals up to this order constrained)
          eta_given : {False, True}, optional
            Set this True if you are passing the marginals in Pr instead of 
            the probabilities
          ic_offset : float, oprtional
            Initial condition offset for the numerical optimisation. If you
            are having trouble getting convergence, try playing with this. 
            Usually making it smaller is effective (ie -0.00001)

        :Returns:
          Psolve : (fdim,)
            probability distribution vector of k-th order maximum entropy
            solution


        """
        if len(Pr.shape) != 1:
            raise ValueError, "Input Pr should be a 1D array"
        if not eta_given and Pr.size != self.fdim:
            raise ValueError, "Input probability vector must have length fdim (m^n)"
        if eta_given:
            if Pr.size != self.dim:
                raise ValueError, "Input eta vector must have length dim (m^n -1)"
        else:
            if Pr.size != self.fdim:
                raise ValueError, "Input probability vector must have length fdim (m^n)"
            if not np.allclose(Pr.sum(), 1.0):
                raise ValueError, "Input probability vector must sum to 1"


        l       = self.order_idx[k].astype(int)
        theta0  = np.zeros(self.order_idx[-1]-self.order_idx[k]-1)
        x0      = np.zeros(l)+ic_offset 
        sf      = self._solvefunc

        jacobian = kwargs.get('jacobian',True)

        Asmall = self.A[:l,:]
        Bsmall = Asmall.T
        if eta_given:
            eta_sampled = Pr[:l]
        else:
            eta_sampled = Asmall*Pr[1:]

        if jacobian:
            self.optout = opt.fsolve(sf, x0, (Asmall,Bsmall,eta_sampled, l), 
                fprime=self._jacobian, col_deriv=1, full_output=1)
        else:
            self.optout = opt.fsolve(sf, x0, (Asmall,Bsmall,eta_sampled, l), 
                full_output=1)

        #self.optout = opt.leastsq(sf, x0, (Asmall,Bsmall,eta_sampled), 
                #full_output=1)
        the_k = self.optout[0]

        print "order: " + str(k) + \
                " ierr: " + str(self.optout[2]) + " - " + self.optout[3]
        print "fval: " + str(np.mean(np.abs(self.optout[1]['fvec']))),
        # extra debug info for jacobian 
        print "nfev: %d" % self.optout[1]['nfev'],
        try:
            print "njev: %d" % self.optout[1]['njev']
        except KeyError:
            print ""
        Psolve = np.zeros(self.fdim)
        Psolve[1:] = self._p_from_theta(np.r_[the_k,theta0])
        Psolve[0] = 1.0 - Psolve.sum()
        return Psolve

    def _solvefunc(self, theta_un, Asmall, Bsmall, eta_sampled, l):
        b = np.exp(Bsmall*theta_un)
        y = eta_sampled - ( (Asmall*b) / (b.sum()+1) )
        return y

    def _jacobian(self, theta, Asmall, Bsmall, eta_sampled, l):
        x = np.exp(Bsmall*theta)
        p = Asmall*x
        q = x.sum() + 1

        J = np.outer(p,p)
        xd = sparse.spdiags(x,0,x.size,x.size,format='csc')
        qdp = (Asmall * xd) * Bsmall
        qdp *= q
        J = J - qdp
        J /= (q*q)

        return J

    def _p_from_theta(self, theta):
        """Internal version - stays in dim space (missing p[0])"""
        pnorm = lambda p: ( p / (p.sum()+1) )
        return pnorm(np.exp(self.A.T*theta))

    def p_from_theta(self, theta):
        """Return full ``fdim`` p-vector from ``fdim-1`` length theta"""
        p = np.zeros(self.fdim)
        p[1:] = self._p_from_theta(theta)
        p[0] = 1.0 - p.sum()
        return p

    def theta_from_p(self, p):
        """Return theta vector from full probaility vector"""
        b = np.log(p[1:]) - np.log(p[0])
        if HAS_UMFPACK:
            # use prefactored matrix
            theta = self.umf.solve(um.UMFPACK_A, self.B, b, autoTranspose=True)
        else:
            theta = spsolve(self.B, b)
        # add theta(0) or not?
        return theta

    def eta_from_p(self, p):
        """Return eta-vector (marginals) from full probability vector"""
        return self.A*p[1:]


def inscol(x,h,n):
    xs = x.shape
    hs = h.shape
 
    if hs[0]==1:    # row vector
        h=h.T
        hs=h.shape

    if n==0:
        y = np.hstack((h,x))
    elif n==xs[1]:
        y = np.hstack((x,h))
    else:
        y = np.hstack((x[:,:n],h,x[:,n:]))

    return y


def order1direct(p,a):
    """Compute first order solution directly for testing"""
    if p.size != a.fdim:
        raise ValueError, "Probability vector doesn't match a.fdim"

    # 1st order marginals
    marg = a.eta_from_p(p)[:a.order_idx[1]]
    # output
    p1 = np.zeros(a.fdim)

    the1pos = lambda x,v: ((v-1)*a.n)+x

    # loop over all probabilities (not p(0))
    for i in range(1,a.fdim):
        Pword = dec2base(np.atleast_2d(i).T,a.m,a.n)

        # loop over each variable
        for j in range(a.n):

            # this value
            x = Pword[0][j]
            if x!=0:
                # this is a normal non-zero marginal
                factor = marg[the1pos(j,x)]
            else:
                # this is a zero-value marginal
                factor = 1 - marg[the1pos(j,np.r_[1:a.m])].sum()

            if p1[i]==0:
                # first entry
                p1[i] = factor
            else:
                p1[i] *= factor

    # normalise
    p1[0] = 1.0 - p1.sum()
    return p1
        
