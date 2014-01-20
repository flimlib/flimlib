/*
 * #%L
 * SLIM-curve package for exponential curve fitting of spectral lifetime data.
 * %%
 * Copyright (C) 2010 - 2014 Gray Institute University of Oxford and Board of
 * Regents of the University of Wisconsin-Madison.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

/* Header file for functions defined in lsqnonneg.c */

/**********************************************************************
   Function GCI_lsqnonneg

   This function solves the non-negative least squares problem:
   minimise |Ax-b| subject to x >= 0 (where |v| is the 2-norm of the
   vector v).

   Arguments: A, an m x n matrix, in the form double A[n][m], so
                 the columns of A are A[0], A[1], ..., A[n-1]
              b, the m-vector

              !!! NB: A and B will both be modified unless preserve is
              !!!     non-zero (see below).

              x, the solution will be placed here
              m ) the dimensions of A
              n )
              preserve, copy A and b before solving the problem
              rnorm, (double *) the value of |Ax-b| with the
                determined x if the function was successful or if the
                iteration count was exceeded.  This can be NULL.
              lambda, an n-vector which will contain the dual vector
                 on completion (that is, the Lagrange multipliers).
                 This can be NULL.

   On exit:   
              The return value will be 0 on success, and negative if
                a problem occurred:
                -1: m > MAX_EQNS or m <= 0
                -2: n > MAX_VARS or n <= 0
                -3: iteration count exceeded: more than 3*n iterations
                    performed
                -4: memory allocation problems

**********************************************************************/

int GCI_lsqnonneg(double **A, double *b, double *x, int m, int n,
		  int preserve, double *rnorm, double *lambda);
