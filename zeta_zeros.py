##### COMPUTATION OF NON-TRIVIAL ZETA ZEROS WITH Re = 1/2 #####
# Determine the imaginary part of the first max_M zeros with Re = 1/2
# and positive imaginary part of the Riemann zeta function.

from time import time
import os
import mpmath as mpm
from multiprocessing import Pool
import parallelization as par

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
    ranges = par.cut_ranges(cores,maxz)
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
    max_M = 20
    n_precision = 2000
    n_cores = 2

    mpm.mp.dps = n_precision

    st = time()
    zz = zeta_zeros(max_M,n_precision,n_cores)
    print(f'Computed zeta zeros after {time()-st} s.')


    if not os.path.exists(f'./M{max_M}_p{n_precision}'):
        os.makedirs(f'./M{max_M}_p{n_precision}')

    with open(f'./M{max_M}_p{n_precision}/NImZetaZero_M{max_M}_p{n_precision}.txt','w') as f:
        for num in zz:
            zstr = mpm.nstr(num,n=n_precision)
            f.write(zstr + '\n')