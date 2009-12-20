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
        int N   # number of subwords
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

def nsb(int Nt, int m, np.ndarray[np.int_t, ndim=1] C,
        verbose=True, var=False):
    """Calculate entropy using C NSB implementation from `Spike Train Analysis
    Toolkit <http://neuroanalysis.org/toolkit/index.html>`_.

    :Parameters:
      Nt : int
        Number of trials.
      m : int
        Dimension of space
      C : (m,) int array
        Vector of counts.
      verbose : {True, False}, optional
        Print warnings from NSB routine.
      var : {False, True}, optional
        Return variance in addition to entropy

    :Returns:
      H : float
        Entropy.
      V : float, optional
        Variance (if requested)

    """
    if C.size != m:
        raise ValueError, "Counts vector must be of length m"
    if C.sum() != Nt:
        raise ValueError, "sum(C) must equal Nt"

    # word list
    cdef np.ndarray[np.int_t, ndim=2] wl 
    wl = np.atleast_2d(np.arange(m, dtype=np.int))
    # word count (cast to double)
    cdef np.ndarray[np.float_t, ndim=1] wc
    wc = C.astype(float)

    # create options structure
    cdef options_entropy opts
    opts.possible_words = -1 
    opts.nsb_precision = 1e-6
    opts.E = 1
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

    # create hist1d structure
    cdef hist1d input
    input.P = Nt
    input.C = m
    input.wordlist = <int **>(wl.data)
    input.wordcnt = <double *>(wc.data)

    entropy_nsb(&input, &opts, entropy)

    print "Entropy: ", entropy[0].value
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

    # free everything
    CFreeEst(entropy, &opts)

    if var:
        return entropy[0].value, entropy[0].ve[0].value
    else:
        return entropy[0].value
    
