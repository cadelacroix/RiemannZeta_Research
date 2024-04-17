##### Computing Zeta Zeros #####

from time import time
import pickle
import mpmath as mpm
import os
from multiprocessing import Pool
import parallelization as par

# vec_zzeros - if rang has one entry, it produces a list with the first 
# rang[0] non-trivial zeros of the Riemann zeta function with a precision 
# of Nprecision. 
# If rang has two entries, it produces the list of zeros between the 
# indices rang[0] and rang[1]

def zeta_zeros_partial(n_precision,*rang):
    mpm.mp.dps = n_precision
    if len(rang) == 1:
        vec = [mpm.zetazero(i).imag for i in range(1,rang[0]+1)]
    elif len(rang) == 2:
        vec = [mpm.zetazero(i).imag for i in range(rang[0]+1,rang[1]+1)]
    return vec


def zeta_zeros(max_M,n_cores):
    ranges=par.cut_ranges(n_cores,max_M)
    tasks_zz = [(n_precision,)+elem for elem in ranges]

    st = time()
    with Pool(n_cores) as pool:
        result = pool.starmap(zeta_zeros_partial,tasks_zz)
    print(f'Computed zeta zeros after {time()-st} s.')

    zz = []
    for item in result:
        zz.extend(item)
    return zz


if __name__ == '__main__':
    max_M = 200
    n_precision = 200
    n_cores = 4

    mpm.mp.dps = n_precision

    zz = zeta_zeros(max_M,n_cores)

    if not os.path.exists(f'./M{max_M}_p{n_precision}'):
        os.makedirs(f'./M{max_M}_p{n_precision}')

    with open(f'./M{max_M}_p{n_precision}/zetazeros_M{max_M}_p{n_precision}_py.txt','wb') as file:
        pickle.dump(zz,file)