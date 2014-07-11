#!/usr/bin/python2.7

###############################################################################
##
##  Original authors: ...
##  Last modified 08/07/2014 (tc)
##
##  Main function for DMRG exercise
##
##  Hamiltonian for the N-site spin-1/2 Heisenberg model with 
##  fixed boundary conditions
##
##  Sum_i=1^N J/2 [S^-_i S^+_i+1 + h.c.] + J^z S^z_i S^z_i+1 - h S^z_i 
##
###############################################################################

import math
import time
import numpy as np
import matplotlib.pyplot as plt

from function1_MPS import *
from function2_MPO import *
from function3_LR import *
from function4_optimization import *

##  Set parameters here
p = {
    'N' : 8  ,   #  number of sites
    'd' : 2,     #  Physical dimension (spin-up or spin-down)
    'D' : 100,   #  Bond dimension
    'J' : 1.,    #  Coupling constant for S^+S^-
    'Jz': 1.,    #  S^z coupling constant
    'h' : 0.,    #  Magnetic field
    'nbr_sweeps' : 4, # Number of sweeps (set to zero to use stopping tolerence instead)
    'pause_frequency' : 25, # number of sweeps between pauses
    'stopping_tolerance': 10e-4  # Specify stopping tolerence
}

if __name__=='__main__':

    print("=============== START RUN ===============")
    print("\tSetting up MPS")
    
    ############
    # note (tc): the following block is not needed in principle, since one
    # could directly use p['N'], but it turned out that in some cases it is
    # safer to do it this way
    nbr_sweeps=p['nbr_sweeps']
    stopping_tolerance=p['stopping_tolerance']
    ############

    ##  INITIALIZE WHOLE THING
    mat_M  = create_random_mps(p['N'], p['d'], p['D'])
    mat_M  = initial_right_normalization(mat_M)
    mat_W  = read_ham_Heis_mpo(p['N'], p['d'], p['J'], p['Jz'], p['h'])
    mat_LR = initial_mat_R(mat_M, mat_W)
    
    ##  KEEP TRACK OF ENERGY
    e_sweeps = []
    e_all = []
    current_energy = 0 # assume real current_energy is never zero
    
    ##  PERFORM SWEEPS
    print("Starting to sweep")
    sweep = 1

    ##   Outer loop over sweeps
    while True:
        old_energy = current_energy

        ##   Right sweep
        print 80 * '-'
        print 'sweep %i - left to right (E=%f)\n' % (sweep, current_energy)
        direction = 'left' #1
        for l in range(1, p['N']):
            current_energy, mat_M = matrix_Heff_optimization(mat_M, mat_W, mat_LR, l)
            e_all.append(current_energy)
            mat_M = normalize_one_site_mps(mat_M, l, direction)
            e_all.append(current_energy)
            mat_LR = update_LR(mat_LR,mat_M,mat_W, l, direction)

        ##   Left sweep
        print 80 * '-'
        print 'sweep %i - right to left (E=%f)\n' % (sweep, current_energy)
        direction = 'right' #-1
        for l in range(p['N'], 1, -1):  
            current_energy, mat_M = matrix_Heff_optimization(mat_M, mat_W, mat_LR, l)
            e_all.append(current_energy)
            mat_M = normalize_one_site_mps(mat_M, l, direction)
            mat_LR = update_LR(mat_LR,mat_M,mat_W, l, direction)

        e_sweeps.append(current_energy)

        ##  MAKE SOME PLOTS
        if sweep % p['pause_frequency'] == 0:
            #print("Plotting current energy profile")
            #plt.plot(e_sweeps)
            #plt.xlabel("SWEEPS")
            #plt.ylabel("ENERGY")
            #plt.title("PROGRESS REPORT")
            #plt.draw()
            #plt.show()
            print
            print 'sweep = %i' % sweep
            print 'current_energy: ', current_energy
            time.sleep(1.0)   #   Show the plot for 3 seconds, then close
            #plt.close()
         
        ## Check for stopping criterion
        
        ## If nbr_sweeps is set to 0 then use only the tolerance stopping criterion
        if p['nbr_sweeps'] == 0 and abs(current_energy - old_energy) < p['stopping_tolerance']:
            ##  YAY! WE HAVE CONVERGED (?)
            print("\nConvergence criterion reached!")
            break;
        elif sweep >= p['nbr_sweeps']:
            print("\nMaximum sweeps reached!")
            break;
        else:
            sweep += 1
            print
    
    ##  FINAL PLOT
    print("Plotting final energy profile")
    plt.plot(e_sweeps)
    plt.xlabel("sweeps")
    plt.ylabel("E")
    plt.title('(N,d,D)=(%i,%i,%i) | (J,Jz,h)=(%.3f,%.3f,%.3f)' % (p['N'], p['d'], p['D'], p['J'], p['Jz'], p['h']))
    plt.savefig("heisenberg_dmrg_energy_sweeps.pdf",bbox_inches='tight')
    plt.close()

    plt.plot(e_all)
    plt.xlabel("steps")
    plt.ylabel("E")
    plt.title('(N,d,D)=(%i,%i,%i) | (J,Jz,h)=(%.3f,%.3f,%.3f)' % (p['N'], p['d'], p['D'], p['J'], p['Jz'], p['h']))
    plt.savefig("heisenberg_dmrg_energy_all.pdf",bbox_inches='tight')
    plt.close()

    print("=============== END RUN ===============")

###############################################################################  
