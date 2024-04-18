##### PARALLELIZATION SUBROUTINES #####

import numpy as np

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


# factor_cores - produces the most "balanced" factorization (f1,f2) of n_cores as
# a product of two integers: n_cores = f1 * f2, with f1 <= f2 and f2/f1 minimal. 
# among all such factorizations.

def factor_cores(n_cores):
    f1 = int(np.floor(n_cores**0.5))
    while n_cores % f1 != 0:
        f1 -= 1
    f2 = int(n_cores/f1)
    return (f1,f2)