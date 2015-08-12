MPS-DMRG-LesHouches2014
=======================

Collaborative work on a simple python code for DMRG on the spin-1/2 Heisenberg chain.

The original code has been written during a tutorial organized by Corinna Kollath, as part of the "4th Les Houches school in computational physics" (http://comp-phys-2014.sciencesconf.org).
Credit is to all the [participants](http://comp-phys-2014.sciencesconf.org/resource/listeparticipants), although I have been modifying the code to some extent, trying to debug.
An original version of the code (withouth my modifications) can be found here: http://www.theory.uni-bonn.de/leshouches-code.

At this stage, the program does not work (the energy does not converge).
Once I get it to a working version, I will include proper credit to all the original authors and write down a minimal documentation.

~~Any help with the debugging is welcome!~~  
More than one year after the school, this collaboration never continued..

### Dimensions
The MPS M is a list of N+1 objects, with the following dimensions (up to truncation when a dimension is larger than D):

- M[0]: dummy object, should never enter the program
- M[1]: 1 x d x d
- M[2]: d x d x (d^2)
- ..
- M[N - 1]:  (d^2) x d x d
- M[N - 1]:  d x d x 1

As for the MPO object W, this is a list of N+1 objects with the following dimensions (no truncation):

- W[0]: 1 x 5 x d x d
- W[1]: 5 x 5 x d x d
- ...
- W[N - 1]: 5 x 5 x d x d
- W[N]: 5 x 1 x d x d

The LR object is again a list of N+1 objects, with the following dimensions (up to truncation):

- LR[0]: 1 x 1 x 1
- LR[1]: d x 5 x d
- LR[2]: (d^2) x 5 x (d^2)
- ...
- LR[N - 2]: (d^2) x 5 x (d^2)
- LR[N - 1]: d x 5 x d
- LR[N]:     1 x 1 x 1

Note that:

- the truncation (given by D) only acts on the dimension on indices `a_(l-1)` and `a_l`, i.e. only for M and LR (not for W);
- M[0] and W[0] are dummy objects that should never be addressed in the program;
- LR[0] and LR[N] are used in the program, and they should be correctly initialized;
- for a site l in the bulk, the corresponding L and R objects are L[l-1] and R[l];
- `dim_b=5` is specific for our MPO (the Heisenberg Hamiltonian for spin 1/2).


### Doubt 1
In the function `initial_right_normalization`, I perform the local right-normalization on all sites starting from the site N (the last one) down to site 2 (the second one). For site 1, I use this specific part of the code:
```python
tot = 0.0
for sigma in range(d):
    Mlast_sigma= numpy.matrix(mat_M[last_site][:, sigma, :])
    tot += (Mlast_sigma * Mlast_sigma.H)
mat_M[last_site] /= numpy.sqrt(tot)
```
I am assuming that this is the correct way to obtain `<psi|psi>=1`, is it true?
