# Zero-based approximation of the Riemann zeta function

This project delves into the *Riemann zeta function*, a cornerstone of number theory due to the remarkable connection between its zeros and the distribution of prime numbers. The exact location of its zeros constitutes the renowned Riemann hypothesis.

This project focuses on a sequence of complex functions conceived by Y. Matiyasevich as efficient approximations of $\zeta$. Our goal is to fully elucidate the nature of the approximation error to uncover new insights into the Riemann zeta function by **analyzing vast datasets**. 

## Overview

The code in this repository –written in Python and Julia– generates the data that defines Matiyasevich's approximations of the Riemann zeta function and then computes the approximation errors. The next section of this README provides a quick introduction to the mathematical objects involved in this code. After that, we give an overview of the various scripts and their key functions.

In order to handle computational instability in our calculations, we need to maintain high precision in our floating-point numbers, generally requiring around $10^5$ or more decimal digits. We achieve this using **arbitrary-precision floating-point arithmetic**, available in Python through the `mpmath` package, and in Julia through its native `BigFloat` type. 

The need for arbitrary precision resulted in significant performance issues, which we addressed by (i) implementing the most computationally intensive parts of the code in **Julia**, and (ii) **parallelizing** whenever feasible. 

The Python function `mpmath.zetazero`, which has no equivalent in Julia, calculates the critical-line zeros of $\zeta$ with a specified precision. While it is possible to call a Python function in Julia using `PyCall`, it does not appear feasible to parallelize a function in this context. This explains our choice of Python for this specific part of the code. 

## A bit of the math behind the code

The Riemann zeta function $\zeta$ is defined as the infinite sum 

$$\zeta(s) = \sum_{n=1}^{+\infty} n^{-s}$$ 

for every complex number $s$ with $\mathrm{Re}(s) > 1$. This function has a pole at $s=1$ and can be extended analytically to the entire complex plane. The function $\zeta$ vanishes at all negative even numbers, known as its *trivial* zeros. The remaining zeros (a.k.a. the *non-trivial* ones) are known to lie within the strip of complex numbers $s$ with $0 < \mathrm{Re}(s) < 1$. The Riemann hypothesis predicts that all the non-trivial zeros lie on the *critical line* $\lbrace s \mid \mathrm{Re}(s) = 1/2 \rbrace$. 

Matiyasevich's approximations of $\zeta$ are finite sums of the form $\sum_{n=1}^N \alpha_n \ n^{-s}$ that *interpolate* $\zeta$ based on a finite number of its critical-line zeros. If we denote by $\rho_1, \rho_2, \rho_3, \ldots$ the critical-line zeros with positive imaginary parts, ordered increasingly by their imaginary parts (see image below), then $\zeta$ vanishes on the complex conjugates $\bar\rho_1, \bar\rho_2, \bar\rho_3, \ldots$ as well. 

<p align="center"><img src="./Images/critical.png" width=300 alt="Critical strip and line with a few non-trivial zeros"></p>

The interpolations are then defined by the formula

$$\Omega_M(s) = \sum_{n=1}^{2M+1} \delta_{M,n} \ n^{-s}$$

where $M$ is a positive integer, $\delta_{M,1} = 1$, and the rest of the coefficients $\delta_{M,n}$ are defined as the unique solutions to the linear system

$$\Omega_M(\rho_1) = \cdots = \Omega_M(\rho_M) = \Omega_M(\bar\rho_1) = \cdots = \Omega_M(\bar\rho_M) = 0.$$

Despite being defined only from a finite number of zeros, Matiyasevich observed that the functions $\Omega_M$ "remember" a great deal of information about the Riemann zeta function and have numerous remarkable properties, as described in the reference [REF]. 

In this project, we focus on the quotient

$$\nu_M(s) := \frac{\Omega_M(s)}{\zeta(s)} = \frac{\sum_{n=1}^{2M+1} \delta_{M,n} \ n^{-s}}{\sum_{n=1}^{+\infty} n^{-s}}.$$

By the Möbius inversion formula, the quotient $\nu_M(s)$ may be written as the infinite sum

$$\nu_M(s) = \sum_{n=1}^{+\infty} \mu_{M,n} \ n^{-s},$$

where the coefficients $\mu_{M,n}$ are determined by the formula

$$\mu_{M,n} = \sum_{\overset{{\displaystyle d=1}}{d \mid n}}^{2M+1} \mu\left(\frac{n}{d}\right) \ \delta_{M,d}.$$

and $\mu$ denotes the Möbius mu function. Finally, we consider the finite truncation

$$\nu_{M,K}(s):= \sum_{n=1}^K \mu_{M,n} \ n^{-s}$$

of $\nu_{M}(s)$. Both $\nu_M(s)$ and $\nu_{M, K}(s)$ quantify the difference between Matiyasevich's interpolations and the Riemann zeta function. 

Through our data-driven approach, we aim to describe the dependencies of the functions $\nu_M(s)$ and $\nu_{M, K}(s)$ on the parameters $M, K, s$, generating a comprehensive mathematical conjecture. 

## Scripts and functions

Before executing any of the scripts described below, parameters should be adjusted in their `main` functions.

### **`zeta_zeros.py`**

This script returns a list of the imaginary parts of the first $M$ critical-line zeros $\rho_1,\rho_2,\ldots,\rho_M$ of the Riemann zeta function (see [image](./Images/critical.png) above), multi-threading the execution of the function `mpmath.zetazero`.

**Parameters**
* `max_M` (int). Maximal index $M$ of the critical-line zeros to be computed.
* `n_dps` (int). Precision in number of digits.
* `n_cores` (int). Number of cores for multi-threading. 

**Output.**
Text file (txt) containing the imaginary parts of the first `max_M` critical-line zeros with precision of `n_dps` number of digits, ordered increasingly and separated by new-line characters `\n`. Its path will be `~/Data/p(n_dps)/ImZetaZero_M(max_M)_p(n_dps).txt`.

**Key functions**

TODO

### **`coef_delta.jl`**

This script computes the coefficients $\delta_{M,n}$ (with $n=1, \ldots, 2M+1$) of Matiyasevich's approximation $\Omega_M$ for all $M = 1$, ..., `max_M`, by solving the linear systems corresponding to all of the approximations simultaneously. It is a multi-thread implementation in Julia of the Gauss algorithm without pivots, based on an algorithm due to Matiyasevich and ___.

**Parameters**

* `max_M` (int). Maximal index of $M$ for the approximation
* `max_computed_M` (int).
* `n_dps` (int). 
* `chunk_size` (int). 

**Input.**
Text file (txt) `~/Data/p(n_dps)/ImZetaZero_M(max_computed_M)_p(n_dps).txt` as output by the script `zeta_zeros.py`.
  
**Output.**
Dictionary whose keys are integers from 1 to `max_M` and the value at key $M$ is a list of the coefficients $\delta_{N,n}$ for $n = 1, ..., 2M+1$ recorded as strings. The dictionary is split as several JSON files of size roughly equal to the parameter `chunk_size`. The path of the JSON files will be `~/Data/p$(n_dps)/CoefDelta_M$(start)-$(end)_p$(n_dps).json` where `start` and `end` determine the range of dictionary keys stored in the file.     

### **`nu_scf.jl`**
