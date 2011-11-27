import scipy
import scipy.sparse as sparse
import numpy as np

# Cython
cimport numpy as np
from libc.math cimport exp
from libc.math cimport sqrt

def _sparse_k(double a2, double la2, double b2, double lb2,
        np.ndarray x1, np.ndarray x2, double eps=3.5e-6, int grad=0,
        double s2=0.0):
    cdef int N = x1.shape[0]
    cdef int M = x2.shape[0]
    cdef double d = 0.0
    cdef double k = 0.0
    cdef double ka = 0.0
    cdef double kb = 0.0
    cdef double a = sqrt(a2)
    cdef double b = sqrt(b2)
    cdef double la = sqrt(la2)
    cdef double s = sqrt(s2)

    K = sparse.lil_matrix((N, M), dtype=np.double)
    if grad:
        gradK = sparse.lil_matrix((N, M, 4), dtype=np.double)

    cdef int idx0 = 0
    cdef int idx1 = 0
    for idx0 in range(0, N):
        if grad:
            gradK[idx0, idx0, 3] = 2*s
        for idx1 in range(0, M):
            d  = x1[idx0] - x2[idx1]
            d  = d*d
            da = d/la2
            db = d/lb2
            ka  = a2 * exp(-0.5*da)
            kb  = b2 * exp(-0.5*db)
            k   = ka + kb
            if k > eps:
                K[idx0, idx1] = k

            if grad:
                k = 2*ka/a
                if k > eps:
                    gradK[idx0,idx1,0] = k
                k = 2*kb/b
                if k > eps:
                    gradK[idx0,idx1,1] = k
                k = ka*d/la/la/la
                if k > eps:
                    gradK[idx0,idx1,2] = k

    if grad:
        return K.tocsc(), gradK.tocsc()
    return K.tocsc()

