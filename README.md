MPS-DMRG-LesHouches2014
=======================

Collaborative work on a simple python code for DMRG on the spin-1/2 Heisenberg chain.

The original code has been written during a tutorial organized by Corinna Kollath, as part of the "4th Les Houches school in computational physics" (http://comp-phys-2014.sciencesconf.org).
Credit is to all the participants, although I have been modifying the code to some extent, trying to debug (but possibly introducing even more bugs).
An original version of the code can be found here: http://www.theory.uni-bonn.de/leshouches-code.

At this stage, the program does not work (the energy does not converge).
Once I get it to a working version, I will include proper credit to all the original authors and write down a minimal documentation.

Any help with the debugging is welcome!

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
