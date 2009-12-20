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

from __future__ import division
import numpy as np
from utils import (prob, _probcount, decimalise, pt_bayescount, 
                   dec2base, ent, malog2)

class BaseSystem:
    """Base functionality for entropy calculations common to all systems"""

    def _calc_ents(self, method, sampling, methods):
        """Main entropy calculation function for non-QE methods"""

        self._sample(method=sampling)
        pt = (method == 'pt') or ('pt' in methods)
        plugin = (method == 'plugin') or ('plugin' in methods)
        nsb = (method == 'nsb') or ('nsb' in methods)
        nsbext = (method == 'nsb-ext') or ('nsb-ext' in methods)
        calc = self.calc

        if (pt or plugin): 
            self._calc_pt_plugin(pt)
        if nsb:
            self._calc_nsb(ext=False)
        if nsbext:
            self._calc_nsb(ext=True)
        if 'HshXY' in calc:
            #TODO: not so efficient since samples PY again
            sh = self._sh_instance()
            sh.calculate_entropies(method=method, 
                                   sampling=sampling, 
                                   methods=methods, calc=['HXY'])
            if pt: 
                self.H_pt['HshXY'] = sh.H_pt['HXY']
            if nsb: 
                self.H_nsb['HshXY'] = sh.H_nsb['HXY']
            if nsbext: 
                self.H_nsbext['HshXY'] = sh.H_nsbext['HXY']
            if plugin or pt: 
                self.H_plugin['HshXY'] = sh.H_plugin['HXY']
        if 'HshX' in calc:
            sh = self._shX_instance()
            sh.calculate_entropies(method=method,
                                   sampling=sampling,
                                   methods=methods, calc=['HX'])
            if pt: 
                self.H_pt['HshX'] = sh.H_pt['HX']
            if nsb: 
                self.H_nsb['HshX'] = sh.H_nsb['HX']
            if nsbext: 
                self.H_nsbext['HshX'] = sh.H_nsbext['HX']
            if plugin or pt: 
                self.H_plugin['HshX'] = sh.H_plugin['HX']
            
        if method == 'plugin':
            self.H = self.H_plugin
        elif method == 'pt':
            self.H = self.H_pt
        elif method == 'nsb':
            self.H = self.H_nsb
        elif method == 'nsb-ext':
            self.H = self.H_nsbext

    def _calc_pt_plugin(self, pt):
        """Calculate direct entropies and apply PT correction if required """
        calc = self.calc
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
            H = -(self.PX*malog2(np.ma.array(self.PiX,copy=False,
                    mask=(self.PiX<=np.finfo(np.float).eps)))).sum(axis=0)
            self.H_plugin['ChiX'] = H
            if pt:
                # no PT correction for ChiX
                self.H_pt['ChiX'] = H
        # for adelman style I(k;spike) (bits/spike)
        if 'HXY1' in calc:
            if self.Y_m != 2:
                raise ValueError, \
                "HXY1 calculation only makes sense for spike data, ie Y_m = 2"
            H = ent(self.PXY[:,1])
            self.H_plugin['HXY1'] = H
            if pt:
                self.H_pt['HXY1'] = H + pt_corr(pt_bayescount(self.PXY[:,1],self.Ny[1])) 
        if 'ChiXY1' in calc:
            if self.Y_m != 2:
                raise ValueError, \
                "ChiXY1 calculation only makes sense for spike data, ie Y_m = 2"
            H = -np.ma.array(self.PXY[:,1]*np.log2(self.PX),copy=False,
                    mask=(self.PX<=np.finfo(np.float).eps)).sum()
            self.H_plugin['ChiXY1'] = H
            if pt:
                # no PT for ChiXY1
                self.H_pt['ChiXY1'] = H
    
    def _calc_nsb(self, ext=False):
        """Calculate NSB corrected entropy
        
        :Parameters:
          ext : {False, True}, optional
            Use external 'nsb-entropy' program vs STAT version

        """
        if ext:
            from utils import nsb_entropy
        else:
            # don't catch exceptions - just fail if can't import
            from statk.wrap import nsb_entropy as _nsb_entropy
            # for debugging
            #nsb_entropy = lambda x,y,z: _nsb_entropy(x,y,z,verbose=True)
            nsb_entropy = _nsb_entropy
        calc = self.calc
        H_nsb = {}
        if 'HX' in calc:
            H = nsb_entropy(self.PX, self.N, self.X_dim) 
            H_nsb['HX'] = H
        if 'HY' in calc:
            H = nsb_entropy(self.PY, self.N, self.Y_dim)
            H_nsb['HY'] = H
        if 'HXY' in calc:
            H = 0.0
            for y in xrange(self.Y_dim):
                H += self.PY[y] * nsb_entropy(self.PXY[:,y], self.Ny[y], self.X_dim) 
            H_nsb['HXY'] = H
        if 'SiHXi' in calc:
            H = 0.0
            for i in xrange(self.X_n):
                H += nsb_entropy(self.PXi[:,i], self.N, self.X_m) 
            H_nsb['SiHXi'] = H
        if 'HiXY' in calc:
            H = 0.0
            for i in xrange(self.X_n):
                for y in xrange(self.Y_dim):
                    H += self.PY[y] * nsb_entropy(self.PXiY[:,i,y], self.Ny[y], self.X_m) 
            H_nsb['HiXY'] = H
        if 'HiX' in calc:
            H = nsb_entropy(self.PiX, self.N, self.X_dim)
            H_nsb['HiX'] = H
        if 'ChiX' in calc:
            print "Warning: No NSB correction applied for ChiX"
            H = -(self.PX*malog2(np.ma.array(self.PiX,copy=False,
                    mask=(self.PiX<=np.finfo(np.float).eps)))).sum(axis=0)
            H_nsb['ChiX'] = H
        if ext:
            self.H_nsbext = H_nsb
        else:
            self.H_nsb = H_nsb

    def calculate_entropies(self, method='plugin', sampling='naive', 
                            calc=['HX','HXY'], **kwargs):
        """Calculate entropies of the system.

        :Parameters:
          method : {'plugin', 'pt', 'qe', 'nsb', 'nsb-ext'}
            Bias correction method to use
          sampling : {'naive', 'kt', 'beta:x'}, optional
            Sampling method to use. 'naive' is the standard histrogram method.
            'beta:x' is for an add-constant beta estimator, with beta value
            following the colon eg 'beta:0.01' [1]_. 'kt' is for the 
            Krichevsky-Trofimov estimator [2]_, which is equivalent to 
            'beta:0.5'.

          calc : list of strs
            List of entropy values to calculate from ('HX', 'HY', 'HXY', 
            'SiHXi', 'HiX', 'HshX', 'HiXY', 'HshXY', 'ChiX', 'HXY1','ChiXY1')

        :Keywords:
          qe_method : {'plugin', 'pt', 'nsb', 'nsb-ext'}, optional
            Method argument to be passed for QE calculation ('pt', 'nsb'). 
            Allows combination of QE with other corrections.
          methods : list of strs, optional
            If present, method argument will be ignored, and all corrections 
            in the list will be calculated. Use to comparing results of 
            different methods with one calculation pass.

        :Returns:
          self.H : dict
            Dictionary of computed values.
          self.H_method : dict
            Dictionary of computed values using 'method'.

        Notes
        -----
        * If the PT method is chosen with outputs 'HiX' or 'ChiX' no bias 
          correction will be performed for these terms.

        References
        ----------
        .. [1] T. Schurmann and P. Grassberger, "Entropy estimation of 
           symbol sequences," Chaos,vol. 6, no. 3, pp. 414--427, 1996.
        .. [2] R. Krichevsky and V. Trofimov, "The performance of universal 
           encoding," IEEE Trans. Information Theory, vol. 27, no. 2, 
           pp. 199--207, Mar. 1981. 


        """
        self.calc = calc
        self.methods = kwargs.get('methods',[])
        for m in (self.methods + [method]):
            if m not in ('plugin','pt','qe','nsb','nsb-ext'):
                raise ValueError, 'Unknown correction method : '+str(m)
        methods = self.methods

        # allocate memory for requested calculations
        if any([c in calc for c in ['HXY','HiXY','HY']]):
            # need Py for any conditional entropies
            self.PY = np.zeros(self.Y_dim)
        if any([c in calc for c in ['HX','HshX','ChiX','ChiXY1']]):
            self.PX = np.zeros(self.X_dim)
        if ('HiX' in calc) or ('ChiX' in calc):
            self.PiX = np.zeros(self.X_dim)
        if any([c in calc for c in ['HXY','HXY1','ChiXY1']]):
            self.PXY = np.zeros((self.X_dim,self.Y_dim))
        if 'SiHXi' in calc:
            self.PXi = np.zeros((self.X_m,self.X_n))
        if ('HiXY' in calc) or ('HiX' in calc):
            self.PXiY = np.zeros((self.X_m,self.X_n,self.Y_dim))
        if 'HshXY' in calc:
            self.Xsh = np.zeros(self.X.shape,dtype=np.int)

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

    def I(self, corr=None):
        """Convenience function to compute mutual information
        
        Must have already computed required entropies ['HX', 'HXY']

        :Parameters:
          corr : str, optional
            If provided use the entropies from this correction rather than
            the default values in self.H
        
        """
        try:
            if corr is not None:
                H = getattr(self,'H_%s'%corr)
            else:
                H = self.H
            I = H['HX'] - H['HXY']
        except (KeyError, AttributeError):
            print "Error: must have computed HX and HXY for" + \
            "mutual information"
            return
        return I

    def Ish(self, corr=None):
        """Convenience function to compute shuffled mutual information
        estimate
        
        Must have already computed required entropies
        ['HX', 'HiXY', 'HshXY', 'HXY']

        :Parameters:
          corr : str, optional
            If provided use the entropies from this correction rather than
            the default values in self.H
        
        """
        try:
            if corr is not None:
                H = getattr(self,'H_%s'%corr)
            else:
                H = self.H
            I = H['HX'] - H['HiXY'] + H['HshXY'] - H['HXY']
        except (KeyError, AttributeError):
            print "Error: must have computed HX, HiXY, HshXY and HXY" + \
                    "for shuffled mutual information estimator"
            return
        return I

    def Ishush(self, corr=None):
        """Convenience function to compute full shuffled mutual information
        estimate
        
        Must have already computed required entropies
        ['HX', 'SiHXi', 'HshX', 'HiXY', 'HshXY', 'HXY']

        :Parameters:
          corr : str, optional
            If provided use the entropies from this correction rather than
            the default values in self.H
        
        """
        try:
            if corr is not None:
                H = getattr(self,'H_%s'%corr)
            else:
                H = self.H
            I = (H['HX'] - H['HshX'] + H['SiHXi'] -
                    H['HiXY'] + H['HshXY'] - H['HXY'])
        except (KeyError, AttributeError):
            print "Error: must have computed HX, HshX, SiHXi, " + \
                "HiXY, HshXY and HXY for shuffled mutual information estimator"
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
        except (KeyError, AttributeError):
            print "Error: must compute SiHXi, HiXY, HiX, ChiX and Ish for Pola breakdown"
        return I

    def Ispike(self):
        """Adelman (2003) style information per spike """
        try:
            I = self.H['ChiXY1'] - self.H['HXY1']
        except (KeyError, AttributeError):
            print "Error: must compute ChiXY1, HXY1 for Ispike"
            return
        return I

    def _qe_ent(self, qe_method, sampling, methods):
        """General Quadratic Extrapolation Function"""
        calc = self.calc
        self._qe_prep()
        N = self.N
        N2 = N/2.0
        N4 = N/4.0

        # full length
        # add on methods to do everything (other than qe) with this one call
        self._calc_ents(qe_method,sampling,methods)
        H1 = np.array([v for k,v in sorted(self.H.iteritems())])
        
        # half length
        H2 = np.zeros(H1.shape)
        half_slices = [(2,0), (2,1)]
        for sl in half_slices:
            sys = self._subsampled_instance(sl)
            sys.calculate_entropies(method=qe_method, sampling=sampling, calc=calc)
            H2 += np.array([v for k,v in sorted(sys.H.iteritems())])
            del sys
        H2 = H2 / 2.0
        
        # quarter length
        H4 = np.zeros(H1.shape)
        quarter_slices = [(4,0), (4,1), (4,2), (4,3)]
        for sl in quarter_slices:
            sys = self._subsampled_instance(sl)
            sys.calculate_entropies(method=qe_method, sampling=sampling, calc=calc)
            H4 += np.array([v for k,v in sorted(sys.H.iteritems())])
            del sys
        H4 = H4 / 4.0

        # interpolation
        Hqe = np.zeros(H1.size)
        for i in xrange(H1.size):
            Hqe[i] = np.polyfit([N4,N2,N],
                        [N4*N4*H4[i], N2*N2*H2[i], N*N*H1[i]], 2)[0]
        keys = [k for k,v in sorted(self.H_plugin.iteritems())]
        self.H_qe = dict(zip(keys, Hqe))


class DiscreteSystem(BaseSystem):
    """Class to hold probabilities and calculate entropies of 
    a discrete stochastic system.

    :Attributes:
      PXY : (X_dim, Y_dim)
        Conditional probability vectors on decimalised space P(X|Y). 
        ``PXY[:,i]`` is X probability distribution conditional on ``Y==i``.
      PX : (X_dim,) 
        Unconditional decimalised X probability.
      PY : (Y_dim,)
        Unconditional decimalised Y probability.
      PXi : (X_m, X_n)
        Unconditional probability distributions for individual X components. 
        ``PXi[i,j] = P(X_i==j)``
      PXiY : (X_m, X_n, Y_dim)
        Conditional probability distributions for individual X compoenents.
        ``PXiY[i,j,k] = P(X_i==j | Y==k)``
      PiX : (X_dim,)
        ``Pind(X) = <Pind(X|y)>_y``

    """

    def __init__(self, X, X_dims, Y, Y_dims, qe_shuffle=True):
        """Check and assign inputs.

        :Parameters:
          X : (X_n, t)  int array
            Array of measured input values. X_n variables in X space, t trials
          X_dims : tuple (n, m)
            Dimension of X (input) space; length n, base m words
          Y : (Y_n, t) int array
            Array of corresponding measured output values. Y_n variables in Y
            space, t trials
          Y_dims : tuple (n ,m)
            Dimension of Y (output) space; length n, base m words
          qe_shuffle : {True, False}, optional
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
        self.N = self.X.shape[1]
        self.Ny = np.zeros(self.Y_dim)
        self.qe_shuffle = qe_shuffle
        self.sampled = False
        self.calc = []
        
    def _sample(self, method='naive'):
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
        if any([c in calc for c in ['HiX','HiXY','HXY','HY']]):
            if self.Y_n > 1:
                d_Y = decimalise(self.Y, self.Y_n, self.Y_m)
            else:
                # make 1D
                d_Y = self.Y.reshape(self.Y.size)

        # unconditional probabilities
        if ('HX' in calc) or ('ChiX' in calc):
            self.PX = prob(d_X, self.X_dim, method=method)
            """test docstring fpr PX"""
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
                            if 'HshXY' in calc:
                                # shuffle
                                #np.random.shuffle(oce)
                                shfoce = np.random.permutation(oce)
                                self.Xsh[j,indx] = shfoce
        # Pind(X) = <Pind(X|Y)>_y
        if ('HiX' in calc) or ('ChiX' in calc):
            # construct joint distribution
            words = dec2base(np.atleast_2d(np.r_[0:self.X_dim]).T,self.X_m,self.X_n)
            PiXY = np.zeros((self.X_dim, self.Y_dim))
            PiXY = self.PXiY[words,np.r_[0:self.X_n]].prod(axis=1)
            # average over Y
            self.PiX = np.dot(PiXY,self.PY)

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

    def _sh_instance(self):
        """Return shuffled instance"""
        # do it like this to allow easy inheritence
        return DiscreteSystem(self.Xsh, self.X_dims, self.Y, self.Y_dims)

    def _shX_instance(self):
        """Return shuffled instance"""
        # do it like this to allow easy inheritence
        # unconditional shuffle
        Xsh_un = np.zeros_like(self.X)
        for i in range(self.X_n):
            shindx = np.random.permutation(self.X.shape[1])
            Xsh_un[i,:] = self.X[i,shindx]
        return DiscreteSystem(Xsh_un, self.X_dims, self.Y, self.Y_dims)

    def _qe_prep(self):
        """QE Preparation"""
        if self.qe_shuffle:
            # need to shuffle to ensure even stimulus distribution for QE
            shuffle = np.random.permutation(self.N)
            # fancy indexing makes a copy
            self.X = self.X[:,shuffle]
            self.Y = self.Y[:,shuffle]

        # ensure trials is a multiple of 4 for easy QE
        rem = np.mod(self.X.shape[1],4)
        if rem != 0:
            self.X = self.X[:,:-rem]
            self.Y = self.Y[:,:-rem]
        self.N = self.X.shape[1]

    def _subsampled_instance(self, sub):
        """Return subsampled instance for QE
        
        sub : tuple (df, i) 
              red - reduction factor (2, 4)
              i - interval

        """
        Nred = self.N / sub[0]
        sl = slice(sub[1]*Nred,(sub[1]+1)*Nred)
        return DiscreteSystem(self.X[:,sl], self.X_dims,
                              self.Y[:,sl], self.Y_dims)


class SortedDiscreteSystem(DiscreteSystem):
    """Class to hold probabilities and calculate entropies of a discrete 
    stochastic system when the inputs are available already sorted 
    by stimulus.

    :Attributes:
      PXY : (X_dim, Y_dim)
        Conditional probability vectors on decimalised space P(X|Y). 
        ``PXY[:,i]`` is X probability distribution conditional on ``Y==i``.
      PX : (X_dim,) 
        Unconditional decimalised X probability.
      PY : (Y_dim,)
        Unconditional decimalised Y probability.
      PXi : (X_m, X_n)
        Unconditional probability distributions for individual X components. 
        ``PXi[i,j] = P(X_i==j)``
      PXiY : (X_m, X_n, Y_dim)
        Conditional probability distributions for individual X compoenents.
        ``PXiY[i,j,k] = P(X_i==j | Y==k)``
      PiX : (X_dim,)
        ``Pind(X) = <Pind(X|y)>_y``

    """
    def __init__(self, X, X_dims, Y_m, Ny):
        """Check and assign inputs.

        :Parameters:
          X : (X_n, t) int array
            Array of measured input values. X_n variables in X space, t trials
          X_dims : tuple (n,m)
            Dimension of X (input) space; length n, base m words
          Y_m : int 
            Finite alphabet size of single variable Y
          Ny : (Y_m,) int array
            Array of number of trials available for each stimulus. This should
            be ordered the same as the order of X w.r.t. stimuli. 
            Y_t.sum() = X.shape[1]

        """
        self.X_dims = X_dims
        self.X_n = X_dims[0]
        self.X_m = X_dims[1]
        self.Y_m = Y_m
        self.X_dim = self.X_m ** self.X_n
        self.Y_dim = self.Y_m 
        self.X = np.atleast_2d(X)
        self.Ny = Ny.astype(float)
        self.N = self.X.shape[1]
        self._check_inputs()
        self.sampled = False
        self.qe_shuffle = True
        self.calc = []

    def _check_inputs(self):
        if (not np.issubdtype(self.X.dtype, np.int)):
            raise ValueError, "Inputs must be of integer type"
        if (self.X.max() >= self.X_m) or (self.X.min() < 0):
            raise ValueError, "X values must be in [0, X_m)"
        if (self.X.shape[0] != self.X_n):
            raise ValueError, "X.shape[0] must equal X_n"
        if (self.Ny.size != self.Y_m):
            raise ValueError, "Ny must contain Y_m elements"
        if (self.Ny.sum() != self.N):
            raise ValueError, "Ny.sum() must equal number of X input trials"

    def _sample(self, method='naive'):
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

        # unconditional probabilities
        if ('HX' in calc) or ('ChiX' in calc):
            self.PX = prob(d_X, self.X_dim, method=method)
        if any([c in calc for c in ['HXY','HiX','HiXY','HY']]):
            self.PY = _probcount(self.Ny,self.N,method)
        if 'SiHXi' in calc:
            for i in xrange(self.X_n):
                self.PXi[:,i] = prob(self.X[i,:], self.X_m, method=method)
            
        # conditional probabilities
        if any([c in calc for c in ['HiXY','HXY','HshXY']]):
            sstart=0
            for i in xrange(self.Y_dim):
                send = sstart+self.Ny[i]
                indx = slice(sstart,send)
                sstart = send
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
                            if 'HshXY' in calc:
                                # shuffle
                                #np.random.shuffle(oce)
                                shfoce = np.random.permutation(oce)
                                self.Xsh[j,indx] = shfoce
        # Pind(X) = <Pind(X|Y)>_y
        if ('HiX' in calc) or ('ChiX' in calc):
            # construct joint distribution
            words = dec2base(np.atleast_2d(np.r_[0:self.X_dim]).T,self.X_m,self.X_n)
            PiXY = np.zeros((self.X_dim, self.Y_dim))
            PiXY = self.PXiY[words,np.r_[0:self.X_n]].prod(axis=1)
            # average over Y
            self.PiX = np.dot(PiXY,self.PY)

        self.sampled = True

    def _sh_instance(self):
        """Return shuffled instance"""
        return SortedDiscreteSystem(self.Xsh, self.X_dims, self.Y_m, self.Ny)

    def _shX_instance(self):
        """Return shuffled instance"""
        # unconditional shuffle
        Xsh_un = np.zeros_like(self.X)
        for i in range(self.X_n):
            shindx = np.random.permutation(self.X.shape[1])
            Xsh_un[i,:] = self.X[i,shindx]
        return SortedDiscreteSystem(Xsh_un, self.X_dims, self.Y_m, self.Ny)

    def _qe_prep(self):
        """QE Preparation"""
        if self.qe_shuffle:
            # need to shuffle to ensure even stimulus distribution for QE
            sstart = 0
            oldX = self.X
            self.X = np.zeros_like(oldX)
            for i in xrange(self.Y_m):
                send = sstart + int(self.Ny[i])
                shuffle = np.random.permutation(int(self.Ny[i]))
                self.X[:,sstart:send] = oldX[:,sstart+shuffle]
                sstart = send
                
    def _subsampled_instance(self, sub):
        """Return subsampled instance for QE
        
        sub : tuple (df, i) 
              red - reduction factor (2, 4)
              i - interval

        """

        # reduce each Y data set 
        slices = []
        Ny_new = np.floor(self.Ny/sub[0]).astype(int)
        sstart = 0
        for i in xrange(self.Y_m):
            send = sstart + int(self.Ny[i])
            sl = slice(sstart + (sub[1] * Ny_new[i]), 
                       sstart + ((sub[1]+1) * Ny_new[i]))
            slices.append(sl)
            sstart = send
        X_new = self.X[:,np.r_[tuple(slices)]]

        return SortedDiscreteSystem(X_new, self.X_dims,
                              self.Y_m, Ny_new)


