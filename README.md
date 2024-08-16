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

$$\nu_M(s) := \frac{\Omega_M(s)}{\zeta(s)} = \sum_{n=1}^{+\infty} \mu_{M,n} \ n^{-s}$$

and its finite truncation

$$\nu_{M,K}(s):= \sum_{n=1}^K \mu_{M,n} \ n^{-s}.$$

Both quantify the difference between Matiyasevich's interpolations and the Riemann zeta function. Through our data-driven approach, we aim to describe the dependencies of the functions $\nu_M(s)$ and $\nu_{M, K}(s)$ on the parameters $M, K, s$, generating a comprehensive mathematical conjecture. 

## Scripts and functions

1_zeta_zeros.py
  Given a positive integer (max_M), this Python script returns a list of the first (max_M) zeros of the Riemann zeta function 
  with real part 1/2 and positive imaginary part. The precision of the zeros (in number of digits) is specified as the variable (n_dps). 
  The variable (n_cores) correspond to the number of cores over which the computation is distributed. 

  The output is a text file that contains the imaginary part of the zeros, ordered increasingly, separated by newline characters ('\n').
  The file will be recorded in the path "~/Data/p(n_dps)/ImZetaZero_M(max_M)_p(n_dps).txt".

  Input parameters (max_M, n_dps, n_cores) must be specified in the "__name__ = '__main__'" section of the code.

  The zeros are computed using the function mpmath.zetazero, built in the mpmath package. 


2_coef_delta.jl
  Assume that the file "~/Data/p(n_dps)/ImZetaZero_M(max_computed_M)_p(n_dps).txt" produced by the script 1_zeta_zeros.py is in place, 
  and let (max_M) be a positive integer less than or equal to (max_computed_M). This Julia script computes for all M = 1, ..., (max_M), 
  the coefficients of the unique finite Dirichlet sum

  \[
    \Omega_M(s) = \sum_{n=1}^N \delta_{N,n} n^{-s}
  \]

  that vanishes on the first M zeros of the Riemann zeta function and on its complex conjugates, with $N = 2M + 1$ and $\delta_{N,1} = 1$.
  The computation is performed with a numerical precision of (n_dps) decimal digits. 

  The output is a dictionary whose keys are integers M = 1, ..., max_M. The value at the key M is a list of length N = 2M + 1 that contains
  the coefficients $\delta_{N,n}$ for n = 1, ..., N. The dictionary is then split in several JSON files of maximal size specified by
  the parameter (chunk_size).
