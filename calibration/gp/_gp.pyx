import scipy
import scipy.sparse as sparse
import numpy as np

# Cython
cimport numpy as np
from libc.math cimport exp

def _sparse_k(double a2, double l2, np.ndarray x1, np.ndarray x2, double chi2max=25.0):
    cdef int N = x1.shape[0]
    cdef int M = x2.shape[0]
    cdef double d = 0.0

    K = sparse.lil_matrix((N, M), dtype=np.double)

    cdef int idx0 = 0
    cdef int idx1 = 0
    for idx0 in range(0, N):
        for idx1 in range(0, M):
            d = x1[idx0] - x2[idx1]
            d *= d/l2
            if d < chi2max:
                K[idx0, idx1] = a2 * exp(-0.5*d)
    return K.tocsc()

