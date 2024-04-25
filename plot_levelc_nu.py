import mpmath as mpm
import numpy as np
from sympy import zeta,N,Rational,log
from sympy.abc import n
from ramanujan.pcf import PCF
import matplotlib.pyplot as plt
import json

mpm.mp.dps = 40_000 # Numerical precision in number of decimal digits
depth = 500

def NN(x):
    return N(x,mpm.mp.dps)

def yuricf(K):
    return PCF(2*K,n*(n+1))

def nu_cf(K):
    return 1-1/(yuricf(K).limit(depth).ratio())

cf_val = []
for i in range(1,401):
    val = nu_cf(int(i))
    cf_val.append(val)
    if i%50 == 0:
        print(f"CF for K={i}: {val}")

f = open('Data/_p10000/nu_levelc_MfloorKdiv2_p10000.json')
data = json.load(f)
nu_val = [mpm.mpmathify(data[f"{i}"][0:1001]) for i in range(2,2600,2)]