import scipy
import scipy.sparse as sparse
import numpy as np

# Cython
cimport numpy as np
from libc.math cimport exp

def _sparse_k(double a2, double la2, double b2, double lb2,
        np.ndarray x1, np.ndarray x2, double eps=3.5e-6):
    cdef int N = x1.shape[0]
    cdef int M = x2.shape[0]
    cdef double d = 0.0
    cdef double k = 0.0

    K = sparse.lil_matrix((N, M), dtype=np.double)

    cdef int idx0 = 0
    cdef int idx1 = 0
    for idx0 in range(0, N):
        for idx1 in range(0, M):
            d  = x1[idx0] - x2[idx1]
            d  = d*d
            da = d/la2
            db = d/lb2
            k  = a2 * exp(-0.5*da) + b2 * exp(-0.5*db)
            if k > eps:
                K[idx0, idx1] = k
    return K.tocsc()

