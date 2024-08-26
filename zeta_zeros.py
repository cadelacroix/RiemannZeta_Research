## COMPUTATION OF CRITICAL-LINE ZEROS ##

import os, numpy as np, mpmath as mpm
from time import time
from multiprocessing import Pool

def zeta_zeros_partial(lo,up,n_precision):
    """
    Computes in a single thread the function mpmath.zetazero for values 
    in range(lo,up) with decimal precision equal to n_precision. 

    Parameters
    ----------
    lo, up, n_precision: int

    Returns
    -------
    list of number type mpmath.mpf
    """
    mpm.mp.dps = n_precision
    return [mpm.zetazero(i).imag for i in range(lo+1,up+1)]


def zeta_zeros(max_M,n_precision=100,n_cores=1):
    """
    Computes zeta_zeros_partial(0,max_M,n_precision) in multithread mode,
    with number of threads equal to n_cores.

    Parameters
    ----------
    max_M, n_precision, n_cores: int

    Returns
    -------
    list of number type mpmath.mpf
    """
    mpm.mp.dps = n_precision
    ranges = np.linspace(0,max_M,n_cores+1,dtype=int)
    tasks = [(lo,up,n_precision) for lo,up in zip(ranges,ranges[1:])]
    with Pool(n_cores) as pool:
        result = pool.starmap(zeta_zeros_partial,tasks)

    return [num for item in result for num in item]


if __name__ == '__main__':
    """
    # If executed, this script exports in a TXT file the result of 
    # zeta_zeros(max_M,n_precision,n_cores), separated by spaces. 
    """
    # Parameters 
    max_M = 100
    n_precision = 1_000
    n_cores = 4

    mpm.mp.dps = n_precision
    st = time()

    # Compute zeta zeros
    zz = zeta_zeros(max_M,n_precision,n_cores)
    print(f'Computed zeta zeros after {time()-st} s.')

    # Check if folder exists
    if not os.path.exists(f'Data/p{n_precision}'):
        os.makedirs(f'Data/p{n_precision}')
    
    # Write file
    filepath = f'Data/p{n_precision}/ImZetaZero_M{max_M}_p{n_precision}.txt'
    with open(filepath,'w') as f:
        for num in zz:
            zstr = mpm.nstr(num,n=n_precision)
            f.write(zstr + '\n')