##### COMPUTATION OF NON-TRIVIAL ZETA ZEROS WITH Re = 1/2 #####
# Determine the imaginary part of the first max_M zeros with Re = 1/2
# and positive imaginary part of the Riemann zeta function.

from time import time
import os, numpy as np, mpmath as mpm
from multiprocessing import Pool


# cut_ranges - subdivides the interval [0,MaxM-1] into n_cores integer 
# subintervals. The output is a list of tuples, each containing the starting
# and the end point of the interval. Used for parallelizing.

def cut_ranges(n_cores,*rang):
    if len(rang) == 1:
        lo = 0; up = rang[0]
    elif len(rang) == 2:
        lo = rang[0]; up = rang[1]
    else:
        raise TypeError('cut_ranges admits one or two integer values for defining the range.')
    
    diff = up - lo
    div = diff // (n_cores) * (n_cores)
    rem = diff % n_cores
    endpts = list(np.linspace(0,div,n_cores+1))
    endpts = [int(item) for item in endpts]
    for ind in range(rem):
        endpts[-1-ind] += rem-ind
    ranges = [(lo+endpts[i],lo+endpts[i+1]) for i in range(len(endpts)-1)]
    return ranges

# zeta_zeros_partial - computation of zeros using mpmath.zetazero. Parameters:
# prec: number of correct decimal digits.
# rang: range of the indexes of zeta zeros to be computed. If it has length one,
# it produces a list of the first rang[0] zeros. If it has length two, it 
# produces the list of zeros between the indices rang[0] and rang[1]. 

def zeta_zeros_partial(prec,*rang):
    mpm.mp.dps = prec
    if len(rang) == 1:
        vec = [mpm.zetazero(i).imag for i in range(1,rang[0]+1)]
    elif len(rang) == 2:
        vec = [mpm.zetazero(i).imag for i in range(rang[0]+1,rang[1]+1)]
    return vec


# zeta_zeros - parallel implementation of zeta_zeros_partial to compute the 
# first non-trivial zeros of zeta. Parameters:
# max_M: number of Riemann zeta zeros to be computed.
# n_precision: number of correct decimal digits.
# n_cores: number of processing cores.

def zeta_zeros(maxz,prec=100,cores=1):
    mpm.mp.dps = prec
    ranges = cut_ranges(cores,maxz)
    tasks_zz = [(prec,)+elem for elem in ranges]
    zz = []

    with Pool(cores) as pool:
        result = pool.starmap(zeta_zeros_partial,tasks_zz)

    for item in result:
        zz.extend(item)
    
    return zz


# If executed, this script exports the imaginary parts of the zeta zeros as 
# a text file, separated by spaces. 

if __name__ == '__main__':
    # Parameters 
    max_M = 1_500
    n_precision = 10_000
    n_cores = 80

    mpm.mp.dps = n_precision

    # Compute zeta zeros
    st = time()
    zz = zeta_zeros(max_M,n_precision,n_cores)
    print(f'Computed zeta zeros after {time()-st} s.')

    # Write zeta zeros
    if not os.path.exists(f'Data/p{n_precision}'):
        os.makedirs(f'Data/p{n_precision}')

    with open(f'Data/p{n_precision}/ImZetaZero_M{max_M}_p{n_precision}.txt','w') as f:
        for num in zz:
            zstr = mpm.nstr(num,n=n_precision)
            f.write(zstr + '\n')