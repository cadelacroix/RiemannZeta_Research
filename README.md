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
