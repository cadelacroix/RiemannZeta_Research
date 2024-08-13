# Zero-based approximation of the Riemann zeta function

This project delves into the *Riemann zeta function*, a cornerstone of number theory due to the remarkable connection between its zeros and the distribution of prime numbers. The exact location of its zeros constitutes the renowned Riemann hypothesis.

This project focuses on a sequence of complex functions conceived by Y. Matiyasevich as efficient approximations of $\zeta$. Our goal is to fully elucidate the nature of the approximation error to uncover new insights into the Riemann zeta function by **analyzing vast datasets**. 

## A bit of math before code

The Riemann zeta function $\zeta$ is defined as the infinite sum 

$$\zeta(s) = \sum_{n=1}^{+\infty} n^{-s}$$ 

for every complex number $s$ with $\mathrm{Re}(s) > 1$. This function has a pole at $s=1$ and can be analytically extended to the rest of the complex plane. The function $\zeta$ vanishes at every negative even number; these are known as its *trivial* zeros. In turn, the non-trivial zeros are known to lie in the so-called *critical strip*, that is, the set of complex numbers $s$ with $0 < \mathrm{Re}(s) < 1$. The Riemann hypothesis predicts that all the non-trivial zeros lie on the *critical line* $\lbrace s \mid \mathrm{Re}(s) = 1/2 \rbrace$. 

![Alt text](https://assets.digitalocean.com/articles/alligator/boo.svg "a title")

Matiyasevich's approximations are finite sums of the form $\sum_{n=1}^N \alpha_n \cdot n^{-s}$ that *interpolate* $\zeta$ based on a finite number of its critical-line zeros. Denoting by $\rho_1, \rho_2, \rho_3, \ldots$ the critical-line zeros with positive imaginary parts, ordered by increasing imaginary parts, 

$$\Omega_M()

## Overview




Computing the Riemann zeta function via finite Dirichlet sum. Based on explorations by Yuri Matiyasevich. 

Scripts:

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
