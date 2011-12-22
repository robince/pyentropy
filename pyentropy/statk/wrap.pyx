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
import numpy as np
cimport numpy as np
cimport cython

# matches #define in toolkit_c.h
DEF MAXSIZE=256

cdef extern from "string.h":
    char *strncpy(char *str1, char *str2, size_t count)
    
cdef extern from "toolkit_c.h":

    struct hist1d:
        int P   # number of trials
        int C   # number of unique words
        int N   # length of word
        int **wordlist # list of words that appear
        double *wordcnt # number of times each word appears

    struct message:
        int i,j,k
        char **status
        char **warnings
        char **errors

    struct nv_pair:
        char name[MAXSIZE]
        double value

    struct estimate:
        char name[MAXSIZE]
        double value
        message *messages
        int E
        nv_pair *extras
        int V
        nv_pair *ve

    struct options_entropy:
        double possible_words
        int possible_words_flag
        double nsb_precision
        int nsb_precision_flag
        int *ent_est_meth
        int E
        int **var_est_meth
        int *V

    estimate *CAllocEst(options_entropy *opts)
    void CFreeEst(estimate *input, options_entropy *opts)
    int entropy_nsb(hist1d *inhist, 
                    options_entropy *opts, estimate *entropy)

def nsb_entropy(np.ndarray[np.float_t, ndim=1] P, int N, int dim, 
        verbose=False, var=False, possible_words='recommended', nsb_precision=1e-6):
    """Calculate entropy using C NSB implementation from `Spike Train Analysis
    Toolkit <http://neuroanalysis.org/toolkit/index.html>`_.

    :Parameters:
      P : (dim,) float array
        Probability vector
      N : int
        Number of trials.
      dim : int
        Dimension of space
      verbose : {False, True}, optional
        Print warnings from NSB routine.
      var : {False, True}, optional
        Return variance in addition to entropy
      possible_words : string, optional or int
        Strategy for choosing total number of possible words
        One of ['recommended', 'unique', 'total', 'possible',
                'min_tot_pos', 'min_lim_tot_pos']
        (default: 'recommended').
        Or an positive integer value, see `here 
        <http://neuroanalysis.org/neuroanalysis/goto.do?page=.repository.toolkit_entropts>`
      nsb_precision : float
        Relative precision for numerical integration (default 1e-6)

    :Returns:
      H : float
        Entropy.
      V : float, optional
        Variance (if requested)

    """
    if P.size != dim:
        raise ValueError, "P vector must be of length dim"
    if np.abs(P.sum()-1.0) > np.finfo(np.float).eps:
        raise ValueError, "sum(P) must equal 1"

    # convert to counts
    cdef np.ndarray[np.float_t, ndim=1] C 
    C = np.zeros(dim, dtype=np.float)
    np.around(P*N, out=C)
    C = C[C>0.5] # only keep observed values

    # create hist1d structure
    cdef hist1d input
    input.P = N # number of trials
    input.C = C.size # number of non-zero counts
    input.N = 1 # ignore any multivariate structure

    # word list
    cdef np.ndarray[np.int_t, ndim=2] wl 
    wl = np.atleast_2d(np.arange(C.size, dtype=np.int))
    input.wordlist = <int **>(wl.data)
    input.wordcnt = <double *>(C.data)

    # create options structure
    cdef options_entropy opts
    # set possible_word strategy
    # put in actual values for 'possible' variants to avoid
    # depending on max_possible_words function
    possible_words_codes = { 'recommended': -1, 'unique': -2, 'total': -3,
                             'possible': dim, 'min_tot_pos': min(N,dim), 
                             'min_lim_tot_pos': min(1e5, N, dim) }
    if not isinstance(possible_words, str):
        # pass option direction if numeric
        opts.possible_words = possible_words
    else:
        # look up code
        opts.possible_words = possible_words_codes[possible_words]

    opts.nsb_precision = nsb_precision
    opts.E = 1 # only 1 entropy method
    cdef int ent_meth, var_meth, nV
    ent_meth = 7 # 'nsb'
    var_meth = 0 # 'nsb-var'
    nV = 1
    cdef int *p_var_meth
    p_var_meth = &var_meth
    opts.ent_est_meth = &ent_meth
    opts.var_est_meth = &p_var_meth
    opts.V = &nV
    
    # create estimate return structure
    cdef estimate *entropy
    entropy = CAllocEst(&opts)
    strncpy(entropy[0].name, "nsb", MAXSIZE)
    if var:
        strncpy(entropy[0].ve[0].name, "nsb_var", MAXSIZE)

    entropy_nsb(&input, &opts, entropy)

    cdef message mess
    mess = entropy[0].messages[0]

    if verbose:
        print "Status"
        print "------"
        for i in range(mess.i):
            print mess.status[i]
        print "Warnings"
        print "--------"
        for j in range(mess.j):
            print mess.warnings[j]
        #print "Extras"
        #print "------"
        #for i in range(entropy.E):
            #print entropy[0].extras[i].name, entropy[0].extras[i].value
    # always print errors
    if mess.k > 0:
        print "Errors"
        print "------"
        for k in range(mess.k):
            print mess.errors[k]

    cdef double H
    cdef double V
    H = entropy[0].value
    V = entropy[0].ve[0].value

    # free everything
    CFreeEst(entropy, &opts)

    if var:
        return H, V
    else:
        return H
