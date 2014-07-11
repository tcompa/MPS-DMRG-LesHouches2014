# last modified:    2014-07-09

import numpy
import scipy.linalg

def matrix_Heff_optimization(mat_M, mat_W, mat_LR, l):

	L = mat_LR[l - 1]
	W = mat_W[l]
	R = mat_LR[l]

        #print '[matrix_Heff_opt] site %i | shapes L/R/W: ' % l, L.shape, R.shape, W.shape

	dim_alm1   = L.shape[0]
	dim_blm1   = L.shape[1]
	dim_bl     = W.shape[1]
	dim_sigmal = W.shape[3]
	dim_al     = R.shape[0]

        # first tensordot
        dim_ind_1  = L.shape[1]
        dim_ind_2  = W.shape[0]
        if dim_ind_1 != dim_ind_2:
            exit('ERROR[matrix_Heff] trying to contract %i-dim and %i-dim indices' % (dim_ind_1, dim_ind_2))
	LW = numpy.tensordot(L, W, axes=([1, 0]))

        # second tensordot
        dim_ind_1  = LW.shape[2]
        dim_ind_2  = R.shape[1]
        if dim_ind_1 != dim_ind_2:
            exit('ERROR[matrix_Heff] trying to contract %i-dim and %i-dim indices' % (dim_ind_1, dim_ind_2))
	LWR = numpy.tensordot(LW, R, axes=([2, 1]))

	Heff = LWR.reshape(dim_sigmal * dim_alm1 * dim_al, dim_sigmal * dim_alm1 * dim_al)

	#Mvec=M.reshape(dim_alm1*dim_sigmal*dim_al,1)
	[eigval, eigvec] = scipy.linalg.eigh(Heff, eigvals=(0, 0))
	mat_M[l] = eigvec.reshape(dim_alm1, dim_sigmal, dim_al) # modified (tc): mat_M --> mat_M[l]
	return eigval[0], mat_M # modified (tc): eigval --> eigval[0]

if __name__=='__main__':

    exit('WARNING: What is the first input of matrix_Heff_optimization?\n         Is it a single tensor or the full list?\nExit')

    dim_sigmal=2
    dim_blm1=5
    dim_bl=5
    dim_alm1=3
    dim_al=3

    L=numpy.arange(dim_alm1*dim_blm1*dim_alm1).reshape(dim_alm1,dim_blm1,dim_alm1)
    W=numpy.arange(dim_blm1*dim_bl*dim_sigmal*dim_sigmal).reshape(dim_blm1,dim_bl,dim_sigmal,dim_sigmal)
    R=numpy.arange(dim_al*dim_bl*dim_al).reshape(dim_al,dim_bl,dim_al)
    mat_M=numpy.arange(dim_alm1*dim_sigmal*dim_al).reshape(dim_alm1,dim_sigmal,dim_al);
    R = 0.0001*R
    L = 0.0001*L
    mat_LR=[L,R]
    l=1
    eigval,mat_M= matrix_Heff_optimization(mat_M,W,mat_LR,1)
    print eigval
    print "The new value of M is"
    print mat_M
