/*
 * #%L
 * FLIMLib package for exponential curve fitting of fluorescence lifetime data.
 * %%
 * Copyright (C) 2010 - 2022 University of Oxford and Board of Regents of the
 * University of Wisconsin-Madison.
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

// Updated file from JG, 2.11.02

/** 
 * FLIMLib - Non-negative Least Squares.

Non-neg least squares library function

   This code is based on code and descriptions in the book
     Lawson, C.L. and Hanson, R.J., Solving Least Squares Problems,
     Prentice-Hall, 1974

   A brief description of the various algorithms used is included
   within this file, but for a full explanation, please refer to the
   above-mentioned book or a similar reference.

   Julian Gilbey, Gray Cancer Institute, September 2002
 *
 * \file GCI_Lsqnonneg.c
 */


/* Include files we need */
#include <stdio.h>   /* for NULL */
#include <stdlib.h>  /* for malloc/free */
#include <string.h>  /* for memcpy */
#include <math.h>    /* for fabs and sqrt */
#include "GCI_Lsqnonneg.h"
#ifndef max
#define max(a,b) ((a)<(b) ? (b) : (a))
#endif

/** The greatest number of variables to be determined.  This is to
   allocate an array in advance to save having to use malloc/free */
#define MAX_VARS 50

/** The greatest number of equations to be used.  Again, this is
   to allocate an array in advance */
#define MAX_EQNS 10000

/** Approximate machine epsilon */
#define EPS 1e-9
/** Tolerance */
#define TOL (100*EPS)

/* Internal function prototypes */
void Householder(int mode, double *v, int p, int l, int m, double *u_p,
				 double *C[], int nC);
void GivensCalc(double v_1, double v_2, double *c, double *s, double *r);
void GivensApply(double c, double s, double *z_1, double *z_2);

/**
   Function Householder.

   This function takes as input a vector, a pivot element and a zero
   point, and applies a Householder transformation to the vector to
   make all of its entries zero from the point specified, pivoting
   around the specified pivot entry, modifying the vector in the
   process, as described below.  This function can also take a
   collection of other vectors and apply this Householder
   transformation to these vectors.

   Here is a brief description of Householder transformations: we are
   looking for an orthogonal transformation Q with

      Qv = ( v_0 v_1 ... v_{p-1} -sigma.sqrt(v_p^2+sum_{i=l}^{m-1} v_i^2)
                 v_{p+1} ... v_{l-1} 0 ... 0)^T =: y

      where sigma = +1 if v_p >= 0 and -1 if v_p < 0.

   Here p is the pivot index and l is the zeroing index, where we
   assume 0 <= p <= l < m (where v is an m-vector, zero-based).  If
   these assumptions are not true, we take Q=I.

   Q is given by  Q = I - 2uu^T/|u|^2  for some vector u, but we do
   not actually need to calculate it.  We do the following:

   Compute: s := -sigma.sqrt(v_p^2+sum_{i=l}^m v_i^2)
            u_i := 0, i=0, ..., p-1
            u_p := v_p - s
            u_i := 0, i=p+1, ..., l-1
            u_i := v_i, i=l, ..., m-1
            b := su_p  (so b = -|u|^2/2)
            Q := ( I_m + b^{-1}uu^T  if b != 0
                 ( I_m               if b = 0

   After this, we place the non-zero components of u and y into the
   storage previously occupied by v; we let v_p := y_p and store u_p
   separately.

   To apply this transformation to the set of m-vectors c_j, getting
   d_j = Qc_j, we do:

            t_j := b^{-1}(u^T c_j)
            d_j := c_j + t_j u

   if b != 0, and d_j := c_j if b = 0.

   This function computes u, b, y:=Qv and d_j := Qc_j for all c_j we
   are given.  If u and y have already been computed, we can apply
   this Q to new c_j vectors.

   Arguments: mode: 1 means calculate u, b, y and then Qc_j
                    2 means u and y are precalculated, just calculate Qc_j
              v, the vector to pivot
              NB this vector is modified during execution of this
                  function, as described above
              p, the pivot
              l, the zeroing point
              m, the length of the vector v
              u_p (double *), the calculated value of u_p
              C, double **C, array of pointers to the c_j vectors
              nC, number of c_j vectors
              If C=NULL or nC=0, the code will not attempt to find any
                  Qc_j.
              NB the calculated d_j := Qc_j vectors will overwrite the
                  C array

   Here is the collected algorithm, as outlined above:

   0  // Start here if mode = 1
   1  s := sqrt(v_p^2 + sum_{i=l}^{m-1} v_i^2)
      (Use the calculation trick described below under Givens
      transformations to avoid calculation errors)
   2  if v_p > 0, s := -s
   3  u_p := v_p - s,  v_p := s
   4  // Start here if mode = 2
   5  b := v_p.u_p
   6  if b=0 or no C, stop
   7  for each c_j do steps 8-10
   8  t_j := (c_j[p].u_p + sum_{i=l}^{m-1} c_j[i].v_i) / b
   9  c_j[p] := c_j[p] + t_j.u_p
   10 for i=l,...,m-1,  c_j[i] := c_j[i] + t_j.v_i

*/
 
void Householder(int mode, double *v, int p, int l, int m, double *u_p,
				 double *C[], int nC)
{
	double vmax, vmax_inv, vtmp, s, b, tj;
	double *cj;
	int i, j;

	if (0>p || p>=l || l>m) return;

	/* find the maximum v_i */
	vmax = fabs(v[p]);
	if (mode == 1) {
		for (i=l; i<m; i++)
			vmax = max(vmax, fabs(v[i]));
		if (vmax == 0) return;  /* oops; not much to do here, really */
		vmax_inv = 1/vmax;  /* only calculate this once; does it matter? */

		/* we now calculate s, which is step 1 in the algorithm */
		vtmp = v[p]*vmax_inv;
		s = vtmp*vtmp;
		for (i=l; i<m; i++) {
			vtmp = v[i]*vmax_inv;
			s += vtmp*vtmp;
		}
		s = vmax*sqrt(s);
		/* this is step 1 done */
		if (v[p]>0) s = -s;
		*u_p = v[p]-s;
		v[p] = s;
	} else {
		if (vmax == 0) return;  /* same test as before, only now we are
								   testing only the v_p element, which is
								   actually y_p */
	}

	/* Now calculate the Qc_j vectors */
	if (nC <= 0 || C==NULL) return;
	b = v[p] * (*u_p);
	if (b >= 0) return;

	for (j=0; j<nC; j++) {
		cj = C[j];
		tj = cj[p] * (*u_p);
		for (i=l; i<m; i++)
			tj += cj[i]*v[i];

		if (tj != 0) {
			tj /= b;
			cj[p] += tj*(*u_p);
			for (i=l; i<m; i++)
				cj[i] += tj*v[i];
		}
	}
}

/**
   GivensCalc.

   The GivensCalc function takes as input a 2-vector, and computes
   the Givens rotation for this vector.

   If the vector is v = (v_1 v_2)^T != 0, then the Givens rotation is
   the 2x2 (special) orthogonal matrix

                G = (  c  s )
                    ( -s  c )

   (with c^2 + s^2 = 1) such that Gv = (r 0)^T, where
   r=|v|=sqrt(v_1^2 + v_2^2).

   This function calculates c, s and r, avoiding overflow and
   underflow issues, and returns these values.  (We can allow r to
   overwrite v_1 or v_2 if we wish.)  The calculation is as follows:

   To compute r=sqrt(x^2 + y^2): 
            t := max {|x|, |y|}
            u := min {|x|, |y|}
            r := ( t.sqrt(1+(u/t)^2)  if t != 0
                 ( 0                  if t = 0

*/

void GivensCalc(double v_1, double v_2, double *c, double *s, double *r)
{
	double w, q;

	if (fabs(v_1) > fabs(v_2)) {
		w = v_2/v_1;
		q = sqrt(1+w*w);
		*c = 1/q;
		if (v_1 < 0) *c = - (*c);
		*s = w*(*c);
		*r = fabs(v_1)*q;
	} else {
		if (v_2 == 0) {
			*c = 1; *s = 0; *r = 0;
		} else {
			w = v_1/v_2;
			q = sqrt(1+w*w);
			*s = 1/q;
			if (v_2 < 0) *s = - (*s);
			*c = w*(*s);
			*r = fabs(v_2)*q;
		}
	}
}

/**
   GivensApply.

   The function GivensApply applies a Givens rotation to a
   2-vector.  No rocket science there.

   Here, then, are the algorithms we use:

   Calculate a Givens rotation for the vector (v_1 v_2)^T
   (GivensCalc):

   Arguments: v_1, v_2, the input vector
              c, s, r, pointers to the answers.  r may point to v_1 or v_2

   1  if |v_1| <= |v_2| go to step 8
   2  w := v_2/v_1
   3  q := sqrt(1+w^2)
   4  c := 1/q
   5  if v_1 < 0, c := -c
   6  s := wc
   7  r := |v_1|.q  and stop
   8  if v_2 != 0, go to step 10
   9  [v_1=v_2=0] c := 1, s:= 0, r := 0, stop
   10 w := v_1/v_2
   11 q := sqrt(1+w^2)
   12 s := 1/q
   13 if v_1 < 0, s := -s
   14 c := ws
   15 r := |v_2|.q  and stop

   Apply the Givens rotation G to the vector (z_1 z_2)^T
   (GivensApply):

   Arguments: c, s, components of the Givens rotation matrix
              z_1, z_2, (double *) the vector to apply it to; these
                will be overwritten!

   1  w := z_1.c + z_2.s
   2  z_2 := -z_1.s + z_2.c
   3  z_1 := w

*/

void GivensApply(double c, double s, double *z_1, double *z_2)
{
	double w;

	w = (*z_1)*c + (*z_2)*s;
	*z_2 = -(*z_1)*s + (*z_2)*c;
	*z_1 = w;
}

#define Amatrix(row,col) A[col][row]

/**
   Function GCI_lsqnonneg.

   This function solves the non-negative least squares problem:
   minimise `|Ax-b|` subject to `x >= 0` (where |v| is the 2-norm of the
   vector v).

   The algorithm is similar to the simplex algorithm, and uses
   Lagrange multipliers.

   1.  Set `P:={}, Z:={1,2,...,n}, x:=0`
   2.  Compute `w := A^T (b-Ax)`
   3.  If `Z={}` or `if w_j <= 0` for all j in Z, stop
   4.  Find t in Z such that `w_t = max { w_j : j in Z }`
   5.  Move index t from Z to P
   6.  Let A_P denote the m x n matrix defined by

                  column j of A_P := ( column j of A   if j in P
                                     ( 0               if j in Z

      Compute the n-vector z as a solution of the regular least squares
      problem minimise `| A_P z - b |`.  Note that only the components
      z_j, j in P, are determined by this problem.  Define z_j:=0 for j in Z.
   7.  If `z_j>0` for all j in P, then set `x:=z` and goto step 2
   8.  Find an index q in P such that
        `x_q/(x_q-z_q) = min { x_j/(x_j-z_j) : z_j <= 0, j in P }`
   9.  Set `alpha := x_q/(x_q-z_q)`
   10. Set `x := x + alpha(z-x)`
   11. Move from set P to set Z all j in P for which `x_j=0`.  Go to step 6.

   On termination, x satisfies x_j>0 for j in P and x_j=0 for j in Z, and
   x is a solution to the least squares problem minimise |A_P z - b|.

   Step 6 (finding the solution to the least squares problem
   `|A_P z - b|`) is performed using a QR factorisation of A_P, which in
   turn is calculated using Householder matrices and Givens rotations,
   and clever tricks for keeping track of the QR factorisation of A_P
   as P changes, as we describe below.
*/

int GCI_lsqnonneg(double **A_orig, double *b_orig, double *x, int m, int n,
				  int preserve, double *rnorm_orig, double *lambda)
{
	double *AA[MAX_VARS], bb[MAX_EQNS], **A, *b;  /* to preserve originals */
	double rnorm2, *rnorm;  /* in case rnorm_orig is NULL */
	double ww[MAX_VARS], z[MAX_EQNS], *w;
	double *C[MAX_VARS];  /* for passing to Householder */
	double wmax, asave, u_p, unorm, sum;  /* miscellany */
	int iter, itmax;
	int index[MAX_VARS], p;
	int i, ii, j, jj, k, l, q, zmax, index_zmax;

	/* Meanings of variables:
	   index[] is an array containing information on the sets Z and P.
       index[0..p-1] is the set P
       index[p..n-1] is the set Z

       z[] is working space
	*/

	/* Check variables */
	if (m > MAX_EQNS || m <= 0) return -1;
	if (n > MAX_VARS || n <= 0) return -2;
	w = (lambda == NULL) ? ww : lambda;
	rnorm = (rnorm_orig == NULL) ? &rnorm2 : rnorm_orig;
	if (preserve) {
		/* We allocate one long array and split it into pieces */
		if ((AA[0] = (double *)malloc((unsigned)(m*n)*sizeof(double))) == NULL)
			return -4;
		for (i=0; i<n; i++) {
			AA[i] = AA[0]+m*i;
			for (j=0; j<m; j++)
				AA[i][j] = A_orig[i][j];
		}
		for (j=0; j<m; j++)
			bb[j] = b_orig[j];
		A = AA;
		b = bb;
	} else {
		A = A_orig;
		b = b_orig;
	}

	iter = 0;
	itmax = 3*n;

	/* Initialise the index and x arrays */
	p = 0;
	for (i=0; i<n; i++) {
		x[i] = 0;
		index[i] = i;
	}

	/* Main loop */
	while (p<n && p<m) {
		/* Compute the dual (negative gradient) vector w=A^T(b-Ax)
		   We are only interested in the components of this vector in Z,
		   so we only calculate these.  It turns out that for these
		   components, we have w_j=(A^T.b)_j, and that (A^T.Ax)_j=0,
		   remembering that A is continually modified.  I really don't
		   yet know why this is true. */

		for (i=p; i<n; i++) {
			j = index[i];
			sum = 0;
			for (k=p; k<m; k++)
				sum += Amatrix(k,j) * b[k];
			w[j] = sum;
		}

		wmax = 0;
		while (wmax == 0) {
			/* We now find the maximum positive w_j. */
			zmax = -1;
			for (i=p; i<n; i++) {
				j = index[i];
				if (w[j]>wmax) {
					wmax = w[j];
					zmax = i;
				}
			}

			if (wmax <= 0) /* all w_j non-positive, so done */
				goto terminate;

			/* We move index[zmax] from Z to P
			   We begin the transformation and check the new diagonal element
			   to avoid near linear dependence */
			index_zmax = index[zmax];
			asave = Amatrix(p,index_zmax);

			Householder(1, A[index_zmax], p, p+1, m, &u_p, NULL, 0);
			/* u has ended up in the index_zmax-th column of A */
			unorm = 0;
			for (k=0; k<p; k++)
				unorm += Amatrix(k,index_zmax)*Amatrix(k,index_zmax);
			unorm = sqrt(unorm);

			/* Is the index_zmax-th column sufficiently independent
			   (whatever that means)? */
			if ((unorm != 0 && fabs(Amatrix(p,index_zmax)/unorm) < TOL) ||
				(unorm == 0 && Amatrix(p,index_zmax) == 0)) {
				/* no, it isn't */
				Amatrix(p,index_zmax) = asave;
				w[index_zmax] = 0;
				wmax = 0;
			} else {
				/* it is: copy b[] into z[], update z[] and solve for ztest,
				   the proposed new value for x[index_zmax] */
				double ztest;
				for (i=0; i<m; i++)
					z[i] = b[i];
				C[0] = z;  /* Need to pass Householder an array of pointers */
				Householder(2, A[index_zmax], p, p+1, m, &u_p, C, 1);
				ztest = z[p] / Amatrix(p,index_zmax);
				if (ztest <= 0) {
					/* nah, this won't do, try a different zmax */
					Amatrix(p,index_zmax) = asave;
					w[index_zmax] = 0;
					wmax = 0;
				}
				/* We've now got an acceptable index to move from Z to P!
				   We'll leave this while (wmax==0) loop now, then. */
			}
		}

		/* And take heed of the fact that we've got our index */
		for (i=0; i<m; i++)
			b[i] = z[i];
		index[zmax] = index[p];
		index[p] = index_zmax;
		p++;

		/* We now apply our Householder transformation to all of the
		   columns of A which correspond to index[p], ..., index[n-1].  We
		   could call the Householder function for each column separately,
		   but it's probably more efficient (fewer function calls) to
		   bundle the columns up into one new temporary array.  Note that
		   we've incremented p! */
		if (p<n) {
			for (i=p; i<n; i++) {
				C[i-p] = A[index[i]];  /* the index[i]-th column of A */
			}
			Householder(2, A[index_zmax], p-1, p, m, &u_p, C, n-p);
		}

		/* ... zero the index_zmax-th column of A below the
		   diagonal, and set w[index_zmax]=0 ... */
		for (i=p; i<m; i++)
			Amatrix(i,index_zmax) = 0;
		w[index_zmax] = 0;

		/* ... and solve the triangular system to solve the least squares
		   problem |A_P.z-b| (step 6 of algorithm), putting the solution
		   into z.  Note also that a copy of b is still stored in z when
		   we begin.

		   We have the matrix A_P which has j-th column being index[j],
		   j=0, ..., p-1.  The j-th column is zero below the diagonal.  So
		   we solve it by back-substitution. */
    
		for (i=p-1; i>=0; i--) {
			j = index[i];
			z[i] /= Amatrix(i,j);
			for (k=0; k<i; k++)
				z[k] -= Amatrix(k,j) * z[i];
		}

		/* This next part is called the "secondary loop" in the book */
		while (1) {
			double alpha, tmp;

			if (++iter > itmax)  /* over the iteration count, oops */
				goto terminate;

			alpha = 2;  /* greater than the values we will be
						   comparing, which will all be <= 1 */
			for (i=0; i<p; i++)
				if (z[i] <= 0) {
					tmp = x[index[i]] / (x[index[i]]-z[i]);  /* This is <= 1 */
					if (tmp<alpha) {
						alpha = tmp;
						q = i;  /* keep track of the minimum index */
					}
				}

			/* Were all z_j>0 for j in P?  If so, set x:=z and break out of
			   this inner loop, returning to the main loop. */
			if (alpha > 1) {
				/* the next four lines of code corresponds to lines 244-247 in
				   the Fortran 66 version */
				for (i=0; i<p; i++) {
					x[index[i]] = z[i];
				}
				break;
			}

			/* Interpolate between old x and new z using 0<alpha<=1 */
			for (i=0; i<p; i++)
				x[index[i]] += alpha * (z[i] - x[index[i]]);

			/* Now we have x[index[q]]=0, so modify A and b and the index
			   arrays to move index[q] and any other i with x[i]=0 from set
			   P to set Z */
      
			while (q >= 0) {
				i = index[q];
				x[i] = 0; /* ensure that rounding errors don't creep in */
				/* now use Givens rotations to fix up A and b */
				for (jj=q+1; jj<p; jj++) {
					double c, s;
					ii = index[jj];
					index[jj-1] = ii;
					GivensCalc(Amatrix(jj-1,ii), Amatrix(jj,ii),
							   &c, &s, &Amatrix(jj-1,ii));
					Amatrix(jj,ii) = 0;

					for (l=0; l<n; l++)
						if (l != ii)
							GivensApply(c, s, &Amatrix(jj-1,l), &Amatrix(jj,l));
					GivensApply(c, s, &b[jj-1], &b[jj]);
				}
				/* We've taken one element out of P */
				p--;
				index[p] = i;

				/* Check that the remaining coefficients in set P are
				   feasible, that is, >=0.  They should be, because of
				   the way alpha was determined.  If they are not,
				   this is due to rounding error, so we set any
				   non-positive elements to zero and move them too
				   from set P to set Z */
				for (q=p-1; q>=0; q--)
					if (x[index[q]] <= 0) break;
				/* now if q>=0, then x[index[q]] <= 0,
				   and if q<0, we are done */
			} /* end of while (q >= 0) loop */

			/* OK, all x[i]=0 have now had their indices moved from P
			   to Z.  We can now solve the least squares problem
			   |Az-b| once again, ready for the next secondary loop.
			   So we simply copy b into z and solve the triangular
			   system, as before.  We could reduce the code size a
			   fraction by having this code at the beginning of the
			   loop, but that would make the loop exit code (checking
			   iter) appear in the middle, which would be a tad
			   confusing. */

			for (i=0; i<m; i++)
				z[i] = b[i];
			for (i=p-1; i>=0; i--) {
				j = index[i];
				z[i] /= Amatrix(i,j);
				for (k=0; k<i; k++)
					z[k] -= Amatrix(k,j) * z[i];
			}
		} /* end of secondary loop */
	} /* end of primary loop */    

 terminate:
	sum = 0;
	if (p < m) {
		for (i=p; i<m; i++)
			sum += b[i] * b[i];
	} else {
		/* this is to get the w (dual) vector correct, in case we ever
		   decide to give it as an output to this function. */
		for (i=0; i<n; i++)
			w[i] = 0;
	}

	(*rnorm) = sqrt(sum);

	if (preserve)
		free(AA[0]);

	return (iter > itmax) ? -3 : 0;
}


// Emacs settings:
// Local variables:
// mode: c
// c-basic-offset: 4
// tab-width: 4
// End:
