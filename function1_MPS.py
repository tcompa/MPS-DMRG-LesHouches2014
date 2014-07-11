# last modified: 2014-07-11

import numpy

def create_random_mps(N, d, D):
    """
    Generates a truncated MPS and fills it with random entries
    Inputs:
      N: number of sites
      d: dimension of the local physical space
      D: maximum bond dimension
    Note:
      the function works for even/odd N, but other parts of the program might require even N
    """
    if N % 2 != 0: exit('ERROR[create_random_mps] this function is general, but keep N even')
    mat_M = [None]
    for l in range(1, N / 2 + 1):
        mat_M.append(numpy.random.uniform(size=(min(D, d ** (l - 1)), d, min(D, d ** l))))
    if N % 2 == 1:
        mat_M.append(numpy.random.uniform(size=(min(D, d ** (N / 2)), d, min(D, d ** (N / 2)))))
    for l in range((N + (N % 2)) / 2 + 1, N + 1):
        mat_M.append(numpy.random.uniform(size=(min(D, d ** (N - l + 1)), d, min(D, d ** (N - l)))))
    return mat_M

def normalize_one_site_mps(mat_M, l, direction):
    if (l == len(mat_M) - 1 and direction == 'left') or (l == 1 and direction == 'right'):
        exit('ERROR[normalize_one_site_mps] (1,right) and (N,left) not allowed (%i, %s)' % (l, direction))
    dim_alm1, d, dim_al = mat_M[l].shape
    if direction == 'right':
        M_l_2d = numpy.reshape(mat_M[l], (dim_alm1, d * dim_al))
        U, S, V = numpy.linalg.svd(M_l_2d, full_matrices=0)
        S = numpy.diag(S)
        mat_M[l] = numpy.reshape(V, (dim_alm1, d, dim_al))
        US = numpy.tensordot(U, S, axes=(1, 0))
        mat_M[l - 1] = numpy.tensordot(mat_M[l - 1], US, (2, 0))
    elif direction == 'left':
        M_l_2d = numpy.reshape(mat_M[l], (dim_alm1 * d, dim_al))
        U, S, V = numpy.linalg.svd(M_l_2d, full_matrices=0)
        S = numpy.diag(S)
        mat_M[l] = numpy.reshape(U, (dim_alm1, d, dim_al))
        SV = numpy.tensordot(S, V, axes=(1, 0))
        mat_M[l + 1] = numpy.tensordot(SV, mat_M[l + 1], (1, 0))
    else:
        exit('ERROR[normalize_one_site_mps] direction must be left/right')
    return mat_M

def initial_right_normalization(mat_M):
    N, d, D = len(mat_M) - 1, mat_M[1].shape[1], mat_M[1].shape[2]
    print '[initial_right_normalization] begin | parameters: N=%i, d=%i, D=%i' % (N, d, D)
    ordered_list_sites = range(N, 1, -1)
    last_site = 1
    # loop over sites N,..,2
    for l in ordered_list_sites:
        print '[initial_right_normalization] right-normalize M[%i]' % l
        mat_M = normalize_one_site_mps(mat_M, l, 'right')
    # last site:
    print '[initial_right_normalization] right-normalize M[%i] (special case)' % last_site
    tot = 0.0
    for sigma in range(d):
        Mlast_sigma= numpy.matrix(mat_M[last_site][:, sigma, :])
        tot += (Mlast_sigma * Mlast_sigma.H)
    mat_M[last_site] /= numpy.sqrt(tot)
    print '[initial_right_normalization] end'
    return mat_M

if __name__=='__main__':
    N = 8
    d = 2
    D = 50
    mat_M = create_random_mps(N, d, D)
    mat_M = initial_right_normalization(mat_M)
