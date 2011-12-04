import scipy
import scipy.sparse as sparse
import numpy as np

# Cython
cimport numpy as np
from libc.math cimport exp
from libc.math cimport sqrt

def _sparse_k(double a2, double la2, double b2, double lb2,
        np.ndarray x1, np.ndarray x2, double eps=3.5e-10, int grad=0,
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
    cdef int Npar = 4
    cdef int idx0 = 0
    cdef int idx1 = 0

    K = sparse.lil_matrix((N, M), dtype=np.double)
    if grad:
        gradK = []
        for idx0 in range(0, Npar):
            gradK.append(sparse.lil_matrix((N, M), dtype=np.double))

    for idx0 in range(0, N):
        if grad:
            gradK[3][idx0, idx0] = 2*s
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
                    gradK[0][idx0,idx1] = k
                k = 2*kb/b
                if k > eps:
                    gradK[1][idx0,idx1] = k
                k = ka*d/la/la/la
                if k > eps:
                    gradK[2][idx0,idx1] = k

    if grad:
        return K.tocsc(), [g.tocsc() for g in gradK]
    return K.tocsc()

