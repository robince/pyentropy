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
import os
from tempfile import NamedTemporaryFile
import subprocess

class DiscreteSystem:
    """Class to hold probabilities and calculate entropies of 
    a discrete stochastic system.

    Attributes
    ----------
    
    PXY[X_dim, Y_dim] :
        Conditional probability vectors on decimalised space P(X|Y). PXY[:,i] 
        is X probability distribution conditional on Y == i.
    PX[X_dim] :
        Unconditional decimalised X probability.
    PY[Y_dim] :
        Unconditional decimalised Y probability.
    PXi[X_m,X_n] :
        Unconditional probability distributions for individual X components. 
        PXi[i,j] = P(X_i==j)
    PXiY[X_m,X_n,Y_dim] :
        Conditional probability distributions for individual X compoenents.
        PXiY[i,j,k] = P(X_i==j | Y==k)
    PiX[X_dim] :
        Pind(X) = <Pind(X|y)>_y

    Methods
    -------

    __init__(X, X_dims, Y, Y_dims) :
        Check and assign inputs.
    calculate_entropies(method='plugin', sampling='naive') :
        Calculate entropies and perform bias correction.
    sample(input, output, method='naive') :
        Sample probabilities of system.

    """

    def __init__(self, X, X_dims, Y, Y_dims, qe_shuffle=True):
        """Check and assign inputs.

        Parameters
        ----------
        X[X_n,t] : int array
            Array of measured input values.
        X_dims : tuple (n,m)
            Dimension of X (input) space; length n, base m words
        Y[Y_n,t] : int array
            Array of corresponding measured output values.
        Y_dims : tutple (n,m)
            Dimension of Y (output) space; length n, base m words
        qe_shuffle : bool (True)
            Set to False if trials already in random order, to skip shuffling
            step in QE. Leave as True if trials have structure (ie one stimuli 
            after another).

        """
        self.X_dims = X_dims
        self.Y_dims = Y_dims
        self.X_n = X_dims[0]
        self.X_m = X_dims[1]
        self.Y_n = Y_dims[0]
        self.Y_m = Y_dims[1]
        self.X_dim = self.X_m ** self.X_n
        self.Y_dim = self.Y_m ** self.Y_n
        self.X = np.atleast_2d(X)
        self.Y = np.atleast_2d(Y)
        self._check_inputs(self.X, self.Y)
        self.N = X.shape[1]
        self.Ny = np.zeros(self.Y_dim)
        self.qe_shuffle = qe_shuffle
        self.sampled = False
        
    def sample(self, method='naive'):
        """Sample probabilities of system.

        Parameters
        ----------

        method : {'naive', 'beta:x', 'kt'}, optional
            Sampling method to use. 'naive' is the standard histrogram method.
            'beta:x' is for an add-constant beta estimator, with beta value
            following the colon eg 'beta:0.01' [1]_. 'kt' is for the 
            Krichevsky-Trofimov estimator [2]_, which is equivalent to 
            'beta:0.5'.

        References
        ----------
        .. [1] T. Schurmann and P. Grassberger, "Entropy estimation of 
           symbol sequences," Chaos,vol. 6, no. 3, pp. 414--427, 1996.
        .. [2] R. Krichevsky and V. Trofimov, "The performance of universal 
           encoding," IEEE Trans. Information Theory, vol. 27, no. 2, 
           pp. 199--207, Mar. 1981. 

        """
        calc = self.calc

        # decimalise
        if any([c in calc for c in ['HXY','HX']]):
            if self.X_n > 1:
                d_X = decimalise(self.X, self.X_n, self.X_m)
            else:
                # make 1D
                d_X = self.X.reshape(self.X.size)
        if any([c in calc for c in ['HiX','HiXY','HXY']]):
            if self.Y_n > 1:
                d_Y = decimalise(self.Y, self.Y_n, self.Y_m)
            else:
                # make 1D
                d_Y = self.Y.reshape(self.Y.size)

        # unconditional probabilities
        if ('HX' in calc) or ('ChiX' in calc):
            self.PX = prob(d_X, self.X_dim, method=method)
        if any([c in calc for c in ['HXY','HiX','HiXY','HY']]):
            self.PY = prob(d_Y, self.Y_dim, method=method)
        if 'SiHXi' in calc:
            for i in xrange(self.X_n):
                self.PXi[:,i] = prob(self.X[i,:], self.X_m, method=method)
            
        # conditional probabilities
        if any([c in calc for c in ['HiXY','HXY','HshXY']]):
            for i in xrange(self.Y_dim):
                indx = np.where(d_Y==i)[0]
                self.Ny[i] = indx.size
                if 'HXY' in calc:
                    # output conditional ensemble
                    oce = d_X[indx]
                    if oce.size == 0:
                        print 'Warning: Null output conditional ensemble for ' + \
                          'output : ' + str(i)
                    else:
                        self.PXY[:,i] = prob(oce, self.X_dim, method=method)
                if any([c in calc for c in ['HiX','HiXY','HshXY']]):
                    for j in xrange(self.X_n):
                        # output conditional ensemble for a single variable
                        oce = self.X[j,indx]
                        if oce.size == 0:
                            print 'Warning: Null independent output conditional ensemble for ' + \
                                'output : ' + str(i) + ', variable : ' + str(j)
                        else:
                            self.PXiY[:,j,i] = prob(oce, self.X_m, method=method)
                            # shuffle
                            np.random.shuffle(oce)
                            self.Xsh[j,indx] = oce
        # Pind(X) = <Pind(X|Y)>_y
        if ('HiX' in calc) or ('ChiX' in calc):
            # average over Y
            PiXi = np.dot(self.PXiY, self.PY)
            # construct joint distribution
            words = dec2base(np.atleast_2d(np.r_[0:self.X_dim]).T,self.X_m,self.X_n)
            self.PiX = PiXi[words,np.r_[0:self.X_n]].prod(axis=1)
                            
        self.sampled = True

    def _check_inputs(self, X, Y):
        if (not np.issubdtype(X.dtype, np.int)) \
        or (not np.issubdtype(Y.dtype, np.int)):
            raise ValueError, "Inputs must be of integer type"
        if (X.max() >= self.X_m) or (X.min() < 0):
            raise ValueError, "X values must be in [0, X_m)"
        if (Y.max() >= self.Y_m) or (Y.min() < 0):
            raise ValueError, "Y values must be in [0, Y_m)"        
        if (X.shape[0] != self.X_n):
            raise ValueError, "X.shape[0] must equal X_n"
        if (Y.shape[0] != self.Y_n):
            raise ValueError, "Y.shape[0] must equal Y_n"
        if (Y.shape[1] != X.shape[1]):
            raise ValueError, "X and Y must contain same number of trials"

    def calculate_entropies(self,  method='plugin', sampling='naive', calc=['HX','HXY'], **kwargs):
        """Calculate entropies of the system.

        Parameters
        ----------
        method : {'plugin', 'pt', 'qe', 'nsb'}
            Bias correction method to use.
        sampling : {'naive', 'kt', 'beta:x'}
            Sampling method to use. See docstring of sampling function for 
            further details
        calc :  {'HX','HY','HXY','SiHXi','HiX','HiXY','HshXY'}
            List of entropies to compute.

        Other Parameters
        ----------------
        qe_method : {'pt', 'nsb'}
            Method argument to be passed for QE calculation. Allows combination
            of QE with other corrections.
        methods : list
            If present, method argument will be ignored, and all corrections in 
            list will be calculated. Used for comparing results of different methods
            with one calculation pass.

        Output
        ------
        self.H : dict
            Dictionary of computed values.

        Notes
        -----
        - If the PT method is chosen with outputs 'HiX' or 'ChiX' no bias 
          correction will be performed for these terms.

        """
        self.calc = calc
        self.methods = kwargs.get('methods',[])
        for m in (self.methods + [method]):
            if m not in ('plugin','pt','qe','nsb'):
                raise ValueError, 'Unknown correction method : '+str(m)
        methods = self.methods

        # allocate memory for requested calculations
        if any([c in calc for c in ['HXY','HiXY','HY']]):
            # need Py for any conditional entropies
            self.PY = np.zeros(self.Y_dim)
        if ('HX' in calc) or ('ChiX' in calc):
            self.PX = np.zeros(self.X_dim)
        if ('HiX' in calc) or ('ChiX' in calc):
            self.PiX = np.zeros(self.X_dim)
        if 'HXY' in calc:
            self.PXY = np.zeros((self.X_dim,self.Y_dim))
        if 'SiHXi' in calc:
            self.PXi = np.zeros((self.X_m,self.X_n))
        if ('HiXY' in calc) or ('HiX' in calc):
            self.PXiY = np.zeros((self.X_m,self.X_n,self.Y_dim))
        if 'HshXY' in calc:
            self.Xsh = np.zeros(self.X.shape,dtype=np.int)

        self.sample(method=sampling)

        if (method == 'qe') or ('qe' in methods):
            # default to plugin method if not specified
            qe_method = kwargs.get('qe_method','plugin')
            if qe_method == 'qe':
                raise ValueError, "Can't use qe for qe_method!"
            self._qe_ent(qe_method,sampling,methods)
            if method == 'qe':
                self.H = self.H_qe
        else:
            self._calc_ents(method, sampling, methods)

    def _calc_ents(self, method, sampling, methods):
        pt = (method == 'pt') or ('pt' in methods)
        plugin = (method == 'plugin') or ('plugin' in methods)
        nsb = (method == 'nsb') or ('nsb' in methods)
        calc = self.calc

        if (pt or plugin): 
            pt_corr = lambda R: (R-1)/(2*self.N*np.log(2))
            self.H_plugin = {}
            if pt: self.H_pt = {}
            # compute basic entropies
            if 'HX' in calc:
                H = ent(self.PX)
                self.H_plugin['HX'] = H
                if pt:
                    self.H_pt['HX'] = H + pt_corr(pt_bayescount(self.PX, self.N))
            if 'HY' in calc:
                H = ent(self.PY)
                self.H_plugin['HY'] = H
                if pt:
                    self.H_pt['HY'] = H + pt_corr(pt_bayescount(self.PY, self.N))
            if 'HXY' in calc:
                H = (self.PY * ent(self.PXY)).sum()
                self.H_plugin['HXY'] = H
                if pt:
                    for y in xrange(self.Y_dim):
                        H += pt_corr(pt_bayescount(self.PXY[:,y], self.Ny[y]))
                    self.H_pt['HXY'] = H
            if 'SiHXi' in calc:
                H = ent(self.PXi).sum()
                self.H_plugin['SiHXi'] = H
                if pt:
                    for x in xrange(self.X_n):
                        H += pt_corr(pt_bayescount(self.PXi[:,x],self.N))
                    self.H_pt['SiHXi'] = H
            if 'HiXY' in calc:
                H = (self.PY * ent(self.PXiY)).sum()
                self.H_plugin['HiXY'] = H
                if pt:
                    for x in xrange(self.X_n):
                        for y in xrange(self.Y_dim):
                            H += pt_corr(pt_bayescount(self.PXiY[:,x,y],self.Ny[y]))
                    self.H_pt['HiXY'] = H
            if 'HiX' in calc:
                H = ent(self.PiX)
                self.H_plugin['HiX'] = H
                if pt:
                    # no PT correction for HiX
                    self.H_pt['HiX'] = H
            if 'ChiX' in calc:
                H = -np.ma.array(self.PX*np.log2(PiX),copy=False,
                        mask=(p<=np.finfo(np.float).eps)).sum(axis=0)
                self.H_plugin['ChiX'] = H
                if pt:
                    # no PT correction for ChiX
                    self.H_pt['ChiX'] = H
        
        if nsb:
            # TODO: 1 external program call if all y have same number of trials
            self.H_nsb = {}
            if 'HX' in calc:
                H = nsb_entropy(self.PX, self.N, self.X_dim)[0] / np.log(2)
                self.H_nsb['HX'] = H
            if 'HY' in calc:
                H = nsb_entropy(self.PY, self.N, self.Y_dim)[0] / np.log(2)
                self.H_nsb['HY'] = H
            if 'HXY' in calc:
                H = 0.0
                for y in xrange(self.Y_dim):
                    H += self.PY[y] * nsb_entropy(self.PXY[:,y], self.Ny[y], self.X_dim)[0] / np.log(2)
                self.H_nsb['HXY'] = H
            if 'SiHXi' in calc:
                # TODO: can easily use 1 call here
                H = 0.0
                for i in xrange(self.X_n):
                    H += nsb_entropy(self.PXi[:,i], self.N, self.X_m)[0] / np.log(2)
                self.H_nsb['SiHXi'] = H
            if 'HiXY' in calc:
                H = 0.0
                for i in xrange(self.X_n):
                    for y in xrange(self.Y_dim):
                        H += self.PY[y] * nsb_entropy(self.PXiY[:,i,y], self.Ny[y], self.X_m)[0] / np.log(2)
                self.H_nsb['HiXY'] = H
            if 'HiX' in calc:
                H = nsb_entropy(self.PiX, self.N, self.X_dim)[0] / np.log(2)
                self.H_nsb['HiX'] = H

        if 'HshXY' in calc:
            #TODO: not so efficient since samples PY again
            sys = DiscreteSystem(self.Xsh, self.X_dims, self.Y, self.Y_dims)
            sys.calculate_entropies(method=method, sampling=sampling, methods=methods, calc=['HXY'])
            if pt: 
                self.H_pt['HshXY'] = sys.H_pt['HXY']
            if nsb: 
                self.H_nsb['HshXY'] = sys.H_nsb['HXY']
            if plugin or pt: 
                self.H_plugin['HshXY'] = sys.H_plugin['HXY']

        if method == 'plugin':
            self.H = self.H_plugin
        elif method == 'pt':
            self.H = self.H_pt
        elif method == 'nsb':
            self.H = self.H_nsb

    def _qe_ent(self, qe_method, sampling, methods):
        calc = self.calc
        rem = np.mod(self.X.shape[1],4)
        if rem != 0:
            self.X = self.X[:,:-rem]
            self.Y = self.Y[:,:-rem]
        N = self.X.shape[1] 
        N2 = N/2.0
        N4 = N/4.0

        if self.qe_shuffle:
            # need to shuffle to ensure even stimulus distribution for QE
            shuffle = np.random.permutation(N)
            self.X = self.X[:,shuffle]
            self.Y = self.Y[:,shuffle]

        # full length
        # add on methods to do everything (other than qe) with this one call
        if 'qe' in methods:
            methods.remove('qe')
        self._calc_ents(qe_method,sampling,methods)
        H1 = np.array([v for k,v in sorted(self.H.iteritems())])
        
        # half length
        sys = DiscreteSystem(self.X[:,:N2], self.X_dims, 
                             self.Y[:,:N2], self.Y_dims)
        sys.calculate_entropies(method=qe_method, sampling=sampling, calc=calc)
        H2a = np.array([v for k,v in sorted(sys.H.iteritems())])
        sys = DiscreteSystem(self.X[:,N2:], self.X_dims, 
                             self.Y[:,N2:], self.Y_dims)
        sys.calculate_entropies(method=qe_method, sampling=sampling, calc=calc)
        H2b = np.array([v for k,v in sorted(sys.H.iteritems())])
        H2 = (H2a + H2b) / 2.0
        
        # quarter length
        sys = DiscreteSystem(self.X[:,:N4], self.X_dims, 
                             self.Y[:,:N4], self.Y_dims)
        sys.calculate_entropies(method=qe_method, sampling=sampling, calc=calc)
        H4a = np.array([v for k,v in sorted(sys.H.iteritems())])
        sys = DiscreteSystem(self.X[:,N4:N2], self.X_dims, 
                             self.Y[:,N4:N2], self.Y_dims)
        sys.calculate_entropies(method=qe_method, sampling=sampling, calc=calc)
        H4b = np.array([v for k,v in sorted(sys.H.iteritems())])
        sys = DiscreteSystem(self.X[:,N2:N2+N4], self.X_dims, 
                             self.Y[:,N2:N2+N4], self.Y_dims)
        sys.calculate_entropies(method=qe_method, sampling=sampling, calc=calc)
        H4c = np.array([v for k,v in sorted(sys.H.iteritems())])
        sys = DiscreteSystem(self.X[:,N2+N4:], self.X_dims, 
                             self.Y[:,N2+N4:], self.Y_dims)
        sys.calculate_entropies(method=qe_method, sampling=sampling, calc=calc)
        H4d = np.array([v for k,v in sorted(sys.H.iteritems())])
        H4 = (H4a + H4b + H4c + H4d) / 4.0
        
        # interpolation
        Hqe = np.zeros(H1.size)
        for i in xrange(H1.size):
            Hqe[i] = np.polyfit([N4,N2,N],
                        [N4*N4*H4[i], N2*N2*H2[i], N*N*H1[i]], 2)[0]
        keys = [k for k,v in sorted(self.H_plugin.iteritems())]
        self.H_qe = dict(zip(keys, Hqe))

    def I(self):
        """Convenience function to compute mutual information"""
        try:
            I = self.H['HX'] - self.H['HXY']
        except KeyError:
            print "Error: must have computed HX and HXY for" + \
            "mutual information"
            return
        return I

    def Ish(self):
        """Convenience function for shuffled mutual information"""
        try:
            I = self.H['HX'] - self.H['HiXY'] + self.H['HshXY'] - self.H['HXY']
        except KeyError:
            print "Error: must have computed HX, HiXY, HshXY and HXY" + \
                    "for shuffled mutual information estimator"
            return
        return I

    def pola_decomp(self):
        """Convenience function for Pola breakdown"""
        I = {}
        try:
            I['lin'] = self.H['SiHXi'] - self.H['HiXY']
            I['sig-sim'] = self.H['HiX'] - self.H['SiHXi']
            I['cor-ind'] = -self.H['HiX'] + self.H['ChiX']
            I['cor-dep'] = self.Ish() - self.H['ChiX'] + self.H['HiXY']
        except KeyError:
            print "Error: must compute SiHXi, HiXY, HiX, ChiX and Ish for Pola breakdown"
        return I
               

class SortedDiscreteSystem(DiscreteSystem):
    """Class to hold probabilities and calculate entropies of 
    a discrete stochastic system when the inputs are available already sorted.

    """

    def __init__(self, input, nt, X_m):
        """Check and assign inputs. 
        
        This input format is compatible with `entropy_tb` from Entropy Toolbox 
        for MATLAB (not released yet).

        Parameters
        ----------
        input[X_n,Ntmax,Y_dim] : int array
            Array of measured input values. For each trial t, and Y value y (Y should be 
            a 1-D space) input[:,t,y] is the vector of X response with values in [0,X_m-1].
        nt[Y_dim] : int array
            Number of trials available for each y.
        X_m : int
            Finite alphabet size of X space.

        """
        self.X_n, self.Ntmax, self.Y_dim = input.shape
        self.X_m = X_m
        self.X_dim = self.X_m ** self.X_n
        self.input = input
        self._check_inputs(input)
        self.N = nt.sum()
        self.Ny = nt
        self.sampled = False

    def sample(self, method='naive'):
        """Sample probabilities of system.

        Parameters
        ----------

        method : {'naive', 'beta:x', 'kt'}, optional
            Sampling method to use. 'naive' is the standard histrogram method.
            'beta:x' is for an add-constant beta estimator, with beta value
            following the colon eg 'beta:0.01' [1]_. 'kt' is for the 
            Krichevsky-Trofimov estimator [2]_, which is equivalent to 
            'beta:0.5'.

        References
        ----------
        .. [1] T. Schurmann and P. Grassberger, "Entropy estimation of 
           symbol sequences," Chaos,vol. 6, no. 3, pp. 414--427, 1996.
        .. [2] R. Krichevsky and V. Trofimov, "The performance of universal 
           encoding," IEEE Trans. Information Theory, vol. 27, no. 2, 
           pp. 199--207, Mar. 1981. 

        """
        calc = self.calc

        # decimalise
        if any([c in calc for c in ['HXY','HX']]):
            d_X = decimalise(self.input, self.X_n, self.X_m).T

        # unconditional probabilities
        if ('HX' in calc) or ('ChiX' in calc):
            self.PX = prob(d_X, self.X_dim, method=method)
        if any([c in calc for c in ['HXY','HiX','HiXY','HY']]):
            self.PY = prob(d_Y, self.Y_dim, method=method)
        if 'SiHXi' in calc:
            for i in xrange(self.X_n):
                self.PXi[:,i] = prob(self.X[i,:], self.X_m, method=method)
            
        # conditional probabilities
        if any([c in calc for c in ['HiXY','HXY','HshXY']]):
            for i in xrange(self.Y_dim):
                indx = np.where(d_Y==i)[0]
                self.Ny[i] = indx.size
                if 'HXY' in calc:
                    # output conditional ensemble
                    oce = d_X[indx]
                    if oce.size == 0:
                        print 'Warning: Null output conditional ensemble for ' + \
                          'output : ' + str(i)
                    else:
                        self.PXY[:,i] = prob(oce, self.X_dim, method=method)
                if any([c in calc for c in ['HiX','HiXY','HshXY']]):
                    for j in xrange(self.X_n):
                        # output conditional ensemble for a single variable
                        oce = self.X[j,indx]
                        if oce.size == 0:
                            print 'Warning: Null independent output conditional ensemble for ' + \
                                'output : ' + str(i) + ', variable : ' + str(j)
                        else:
                            self.PXiY[:,j,i] = prob(oce, self.X_m, method=method)
                            # shuffle
                            np.random.shuffle(oce)
                            self.Xsh[j,indx] = oce
        # Pind(X) = <Pind(X|Y)>_y
        if ('HiX' in calc) or ('ChiX' in calc):
            # average over Y
            PiXi = np.dot(self.PXiY, self.PY)
            # construct joint distribution
            words = dec2base(np.atleast_2d(np.r_[0:self.X_dim]).T,self.X_m,self.X_n)
            self.PiX = PiXi[words,np.r_[0:self.X_n]].prod(axis=1)
                            
        self.sampled = True


    


def prob(x, n, method='naive'):
    """Sample probability of integer sequence.

    Parameters
    ----------
    x : int array
        integer input sequence
    n : int
        dimension of input sequence (max(x)<r)
    method: {'naive', 'kt', 'beta:x'}
        Sampling method to use. 

    Returns
    -------
    Pr : float array
        array representing probability distribution
        Pr[i] = P(x=i)

    """
    if (not np.issubdtype(x.dtype, np.int)): 
        raise ValueError, "Input must be of integer type"

    P = np.bincount(x).astype(np.float)
    r = P.size
    if r < n:   # resize if any responses missed
        P.resize((n,))
        P[n:]=0

    if method.lower() == 'naive':
        # normal estimate
        P /= x.size
    elif method.lower() == 'kt':
        # KT (constant addition) estimate
        P = (P + 0.5) / (x.size + (n/2.0))
    elif method.split(':')[0].lower() == 'beta':
        beta = float(method.split(':')[1])
        # general add-constant beta estimate
        P = (P + beta) / (x.size + (beta*r))
    else:
        raise ValueError, 'Unknown sampling method: '+str(est)

    return P


def decimalise(x, n, m):
    """Decimalise discrete response.

    Parameters
    ----------
    x[n,t]: int array
        Vector of samples. Each sample t, is a length-n base-m word.
    n, m : int
        Dimensions of space.

    """
    #TODO: Error checking?
    powers = m**np.arange(0,n,dtype=int)
    #self.powers = self.powers[::-1]
    d_x = np.dot(x.T,powers).astype(int)

    return d_x


ent = lambda p: -np.ma.array(p*np.log2(p),copy=False,
            mask=(p<=np.finfo(np.float).eps)).sum(axis=0)


def pt_bayescount(Pr, Nt):
    """Compute the support for analytic bias correction
    
    Pr - probability
    Nt - number of trials
    
    """
    
    eps = np.finfo(np.float).eps
    # dimension of space
    dim = Pr.size

    Rs = (Pr>eps).sum()

    if Rs < dim:
        Rs_x = Rs - np.exp( Nt*np.log(1.0 - Pr + eps) )[ (Pr>eps) & (Pr<1) ].sum()
        delta_N_prev = dim
        delta_N = np.abs(Rs - Rs_x)
        xtr = 0
        while (delta_N < delta_N_prev) and ((Rs+xtr)<Nt):
            xtr = xtr+1
            Rs_x = 0.0
            gg = xtr*(1.0 - ( (Nt/(Nt+Rs))**(1.0/Nt)))
            qc_x = (1-gg) * (Pr*Nt+1) / (Nt+Rs)
            Rs_x = (1.0 - np.exp(Nt*np.log(1.0-qc_x)))[Pr>eps].sum()
            qc_x = gg / xtr
            Rs_x = Rs_x + xtr*(1.0 - np.exp( Nt * np.log(1.0 - qc_x)))
            delta_N_prev = delta_N
            delta_N = np.abs(Rs - Rs_x)
        Rs = Rs + xtr - 1
        if delta_N < delta_N_prev:
            Rs = Rs + 1

    return Rs


def nsb_entropy(P, N, dim):
    """Calculate NSB entropy of a probability distribution using
    external nsb-entropy program.

    Inputs:
    P - probability distribution vector
    N - total number of trials
    dim - full dimension of space
    
    """
    
    freqs = np.round(P*N)
    tf = NamedTemporaryFile(mode='w',suffix='.txt')

    # write file header
    tf.file.write("# type: scalar\n")
    tf.file.write(str(dim) + "\n")
    tf.file.write("# rows: 1\n")
    tf.file.write("# columns: " + str(freqs.sum().astype(int)) + "\n")

    # write data
    for i in xrange(freqs.size):
        tf.file.write(freqs[i]*(str(i)+" "))
    tf.file.write("\n")
    tf.file.close()

    # run nsb-entropy application
    subprocess.call(["nsb-entropy","-dpar","-iuni","-cY","-s1","-e1",tf.name[:-4]], 
                    stdout=open(os.devnull),stderr=open(os.devnull))

    # read results
    dir, fname = os.path.split(tf.name)
    out_fname = os.path.splitext(fname)[0]+"_uni_num"+str(dim)+"_mf1f0_1_entr.txt"
    out_fname = os.path.join(dir,out_fname)
    fd = open(out_fname,mode='r')
    results = fd.readlines()
    fd.close()
    os.remove(out_fname)

    H = float(results[15].split(' ')[0])
    dH = float(results[20].split(' ')[0])

    return [H, dH]


def dec2base(x, b, digits):
    """Convert decimal value to a row of values representing it in a 
    given base."""

    xs = x.shape
    if xs[1] != 1:
        raise ValueError, "Input x must be a column vector!"

    power = np.ones((xs[0],1)) * (b ** np.c_[digits-1:-0.5:-1,].T)
    x = np.tile(x,(1,digits))

    y = np.floor( np.remainder(x, b*power) / power )

    return y.astype(int)

