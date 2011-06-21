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
"""Utility functions for working with discrete probability
distributions. These functions are exposed in the top-level pyentropy
namespace.

"""
from __future__ import division
import numpy as np
from tempfile import NamedTemporaryFile
import os
import subprocess
from numpy.ma.core import _MaskedUnaryOperation, _DomainGreater
import numpy.core.umath as umath

malog2 = _MaskedUnaryOperation(umath.log2, 1.0, _DomainGreater(0.0))

def ent(p):
    mp = np.ma.array(p,copy=False,mask=(p<=np.finfo(np.float).eps))
    return -(mp*malog2(mp)).sum(axis=0)


def prob(x, m, method='naive'):
    """Sample probability of integer sequence.

    :Parameters:
      x : int array
        integer input sequence
      m : int
        alphabet size of input sequence (max(x)<m)
      method: {'naive', 'kt', 'beta:x','shrink'}
        Sampling method to use. 

    :Returns:
      Pr : float array
        array representing probability distribution
        Pr[i] = P(x=i)

    """
    if not np.issubdtype(x.dtype, np.integer): 
        raise ValueError, "Input must be of integer type"
    if x.max() > m-1:
        raise ValueError, "Input contains values that are too large"
    
    C = np.bincount(x)
    if C.size < m:   # resize if any responses missed
        C.resize((m,))
    return _probcount(C, x.size, method)


def _probcount(C, N, method='naive'):
    """Estimate probability from a vector of bin counts
    
    :Parameters:
      C : int array
        integer vector of bin counts
      N : int
        number of trials
      method: {'naive', 'kt', 'beta:x','shrink'}
        Sampling method to use. 

    """
    N = float(N)
    if method.lower() == 'naive':
        # normal estimate
        P = C/N
    elif method.lower() == 'kt':
        # KT (constant addition) estimate
        P = (C + 0.5) / (N + (C.size/2.0))
    elif method.lower() == 'shrink':
        # James-Stein shrinkage
        # http://www.strimmerlab.org/software/entropy/index.html
        Pnaive = C/N
        target = 1./C.size
        lam = _get_lambda_shrink(N, Pnaive, target)
        P = (lam * target) + ((1 - lam) * Pnaive)
    elif method.split(':')[0].lower() == 'beta':
        beta = float(method.split(':')[1])
        # general add-constant beta estimate
        P = (C + beta) / (N + (beta*C.size))
    else:
        raise ValueError, 'Unknown sampling method: '+str(est)
    return P


def _get_lambda_shrink(N, u, target):
    """Lambda shrinkage estimator"""
    # *unbiased* estimator of variance of u
    varu = u*(1-u)/(N-1)
    # misspecification
    msp = ((u-target)**2).sum()

    # estimate shrinkage intensity
    if msp == 0:
        lam = 1.
    else:
        lam = (varu/msp).sum()
        
    # truncate
    if lam > 1:
        lam = 1 
    elif lam < 0:
        lam = 0

    return lam


def pt_bayescount(Pr, Nt):
    """Compute the support for analytic bias correction using the 
    Bayesian approach of Panzeri and Treves (1996)
    
    :Parameters:
      Pr : 1D aray
        Probability vector
      Nt : int
        Number of trials

    :Returns:
      R : int
        Bayesian estimate of support
    
    """
    
    # dimension of space
    dim = Pr.size

    # non zero probs only
    PrNZ = Pr[Pr>np.finfo(np.float).eps]
    Rnaive = PrNZ.size
    
    R = Rnaive
    if Rnaive < dim:
        Rexpected = Rnaive - ((1.0-PrNZ)**Nt).sum()
        deltaR_prev = dim
        deltaR = np.abs(Rnaive - Rexpected)
        xtr = 0.0
        while (deltaR < deltaR_prev) and ((Rnaive+xtr)<dim):
            xtr = xtr+1.0
            Rexpected = 0.0
            # occupied bins
            gamma = xtr*(1.0 - ((Nt/(Nt+Rnaive))**(1.0/Nt)))
            Pbayes = ((1.0-gamma) / (Nt+Rnaive)) * (PrNZ*Nt+1.0)
            Rexpected = (1.0 - (1.0-Pbayes)**Nt).sum()
            # non-occupied bins
            Pbayes = gamma / xtr
            Rexpected = Rexpected + xtr*(1.0 - (1.0 - Pbayes)**Nt)
            deltaR_prev = deltaR
            deltaR = np.abs(Rnaive - Rexpected)
        Rnaive = Rnaive + xtr - 1.0
        if deltaR < deltaR_prev:
            Rnaive += 1.0
    return Rnaive


def nsb_entropy(P, N, dim):
    """Calculate NSB entropy of a probability distribution using
    external nsb-entropy program.

    Required `nsb-entropy` installed on system path.

    :Parameters:
      P : 1D array
        Probability distribution vector
      N : int
        Total number of trials
      dim : int 
        Full dimension of space
    
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
    given base.
    
    :Parameters:
      x : (t,) or (t,1) int array
        Array of decimilised values
      b : int
        Base for convergence (finite alphabet size)
      digits : int
        Length of output word for each trial

    :Returns:
     y : (t, digits)

    """
    if not np.issubdtype(x.dtype, np.integer):
        raise ValueError, "Input x must be integer"
    if len(x.shape)==1:
        # 1D vector
        x = np.reshape(x,(x.size,1))
    xs = x.shape
    if xs[1] != 1:
        raise ValueError, "Input x must be a 1D array or column vector!"

    power = np.ones((xs[0],1)) * (b ** np.c_[digits-1:-0.5:-1,].T)
    x = np.tile(x,(1,digits))
    y = np.floor( np.remainder(x, b*power) / power )
    return y.astype(int)

def base2dec(x, b):
    """Convert base-b words to decimal values.

    :Parameters:
      x : (t, n) int array
        Array of t length-n base-b words
      b : int
        Base (size of finite alphabet)

    :Returns:
      d_x: (t,)
        Array of decimalised values
    
    Note, this is the same as decimalise except input x is ordered 
    differently (here x[t,n] - ie columns are trials).
    
    """
    xs = x.shape
    z = b**np.arange((xs[1]-1),-0.5,-1)
    d_x = np.dot(x, z)
    return d_x.astype(int)


def decimalise(x, n, b):
    """Convert base-b words to decimal values

    :Parameters:
      x : (n, t) int array
        Array of t length-n base-b words
      b : int
        Base (size of finite alphabet)

    :Returns:
      d_x: (t,)
        Array of decimalised values

    """
    if x.shape[0] != n or x.max() > b-1:
        raise ValueError, "Input vector x doesnt match parameters"
    powers = b**np.arange(n-1,-0.5,-1)
    d_x = np.dot(x.T,powers).astype(int)
    return d_x


def quantise(input, m, uniform='sampling', minmax=None,
             centers=True):
    """ Quantise 1D input vector into m levels (unsigned)

    :Parameters:
      uniform : {'sampling','bins'}
        Determine whether quantisation is uniform for sampling (equally 
        occupied bins) or the bins have uniform widths
      minmax : tuple (min,max)
        Specify the range for uniform='bins' quantisation, rather than using
        min/max of input
      centers : {True, False}
        Return vector of bin centers instead of bin bounds

    """
    bin_centers = np.zeros(m)
    if uniform == 'sampling':
        #bin_numel = np.round(input.size/m) - 1
        bin_numel = np.floor(input.size/m)
        r = input.size - (bin_numel*m) 

        stemp = input.copy()
        stemp.sort(axis=0)
        # original method
        #bin_bounds = stemp[bin_numel:-bin_numel+1:bin_numel]
        # more uniform method
        idx = np.arange(bin_numel, bin_numel*m, bin_numel, dtype=np.int)
        idx[0:r] = idx[0:r] + np.arange(1,r+1,dtype=np.int)
        idx[r:] = idx[r:] + r
        bin_bounds = stemp[idx]

        if centers:
            # calculate center for each bin
            bin_centers[0] =  (bin_bounds[0]+stemp[0]) / 2.0        
            for i in range(1,m-1):
                bin_centers[i] = (bin_bounds[i]+bin_bounds[i-1])/2.0
            bin_centers[m-1] = (stemp[-1]+bin_bounds[-1]) / 2.0
    elif uniform == 'bins':
        if minmax is not None:
            min, max = minmax
        else:
            min, max = input.min(), input.max()
        drange = float(max) - float(min)
        bin_width = drange / float(m)
        bin_bounds = np.arange(1,m,dtype=float)
        bin_bounds *= bin_width
        bin_bounds += min
        if centers:
            bin_centers = r_[bin_bounds - (bin_width/2.0), bin_bounds[-1]+(bin_width/2.0)]
    else:
        raise ValueError, "Unknown value of 'uniform'"

    q_value = np.digitize(input, bin_bounds)

    if centers:
        # bin centers
        return q_value, bin_bounds, bin_centers
    else:
        return q_value, bin_bounds

def quantise_discrete(input, m):
    """Re-bin an already discretised sequence (eg of integer counts)

    Input should already be non-negative integers

    """
    # astype forces a copy even if already int
    X = input.astype(np.int)
    if (X.min() < 0) or not np.all(np.asarray(X,dtype=np.int)==X):
        raise ValueError, "Expecting non-negative integer input"

    if input.max() < m:
        # nothing to do
        return input

    # rebinning algorithm
    # get bincounts now to determine smallest bins
    # merge smallest bin to smallest neighbouring bin
    # (maintain continuity) until we have the right number
    counts = list(np.bincount(X))
    Nbins = len(counts)
    labels = list(np.r_[0:Nbins])

    def merge_bins(a,b):
        """Merge bin a into bin b"""
        counts[b] = counts[b] + counts[a]
        X[X==labels[a]] = labels[b]
        counts.pop(a)
        labels.pop(a)

    while Nbins > m:
        cidx = np.argsort(counts)
        # smallest one
        si = cidx[0]
        # if its at the edges can only merge one way
        if si == 0:
            merge_bins(si,1)
        elif si == len(counts)-1:
            merge_bins(si, si-1)
        else:
            # merge to the smallest neighbour
            target = [si-1, si+1][np.argmin([counts[si-1], counts[si+1]])]
            merge_bins(si, target)
        Nbins = Nbins - 1

    # relabel
    newlabels = range(len(labels))
    for i in range(len(labels)):
        if newlabels[i] != labels[i]:
            # only reassign if necessary
            X[X==labels[i]] = newlabels[i]

    return X


