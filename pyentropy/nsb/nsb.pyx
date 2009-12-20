import numpy as np
cimport numpy as np
cimport cython

cdef extern from "toolkit_c.h":

    struct hist1d:
        int P
        int C
        int N
        int** wordlist
        double* wordcount

    #struct estimate:

    #struct message:

    #struct nv_pair:

    #struct options_entropy:

def test(a,l):
    cdef hist1d x
    x.P = a
    print x.P
        

