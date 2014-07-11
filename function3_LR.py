# last modified: 2014-07-10
# description:   initializes and updates the LR objects

import numpy

def update_LR(mat_LR, mat_M, mat_W, l, direction):
    if l == 0:
        exit('ERROR[update_LR]: l=0 is not accepted as an argument (direction=%s)' % direction)
    if (l == (len(mat_M) - 1) and direction == 'left') or (l == 1 and direction == 'right'):
        exit('ERROR[update_LR]: (1,left) and (N,right) forbidden')

    # construct the MWM object, which has free indices (alm1, al, blm1, bl, alm1p, alp)
    # note: this is local, it does not depend on whether we go left or right
    M = mat_M[l]
    W = mat_W[l]
    # contract sigma_l and sigma_l_prime, to construct MWM
    MW  = numpy.tensordot(M,  W, axes=(1, 2))  # contract sigma_l
    MWM = numpy.tensordot(MW, M, axes=(4, 1)) # contract sigma_l_prime
    # do the last contraction, depending on the direction
    if direction == 'right':
        mat_LR[l - 1] = numpy.tensordot(MWM, mat_LR[l], ((1, 3, 5), (0, 1, 2)))
    elif direction == 'left':
        mat_LR[l] = numpy.tensordot(MWM, mat_LR[l - 1], ((0, 2, 4), (0, 1, 2)))
    else:
        exit('ERROR [update_LR]: \'direction\' has to be right/left')
    return mat_LR


def initial_mat_R(mat_M, mat_W):
    N, d, dim_b = len(mat_W) - 1, mat_W[1].shape[2], mat_W[1].shape[1]
    print '[initial_mat_R] begin (N=%i, d=%i, dim_b=%i)' % (N, d, dim_b)
    if N % 2 != 0: exit('ERROR[initial_mat_R] not yet generalized to odd N')

    # define list structure for mat_LR
    mat_LR = [numpy.ones(shape=(1, 1, 1))]
    for l in range(1, N):
        mat_LR.append(None)
    mat_LR.append(numpy.ones(shape=(1, 1, 1)))

    # construct all R's
    for l in range(N, 1, -1):
        mat_LR = update_LR(mat_LR, mat_M, mat_W, l, 'right')
        print '[update_LR] l=%i, direction=right' % l
    print '[initial_mat_R] end'

    return mat_LR

if __name__=='__main__':

    from function1_MPS import *
    N = 8
    D = 10
    d = 2
    mat_M  = create_random_mps(N, d, D)

    dim_b = 5
    mat_W = [None]
    mat_W.append(numpy.zeros(shape=(1, dim_b, d, d)))
    for site in range(1, N - 1):
        mat_W.append(numpy.zeros(shape=(dim_b, dim_b, d, d)))
    mat_W.append(numpy.zeros(shape=(dim_b, 1, d, d)))
    
    initial_mat_R(mat_M, mat_W)
