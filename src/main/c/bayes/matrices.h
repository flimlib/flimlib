/*-
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
/* mrowley - creation 070904 */


#define USE_COOLEN_VECTOR


#ifdef USE_COOLEN_VECTOR

 float *Bayes_vector(int nl,int nh);
 void free_Bayes_vector(float *v,int nl,int nh);
 double *Bayes_dvector(int nl,int nh);
 void free_Bayes_dvector(double *v,int nl,int nh);
 int *Bayes_ivector(int nl,int nh);
 void free_Bayes_ivector(int *v,int nl,int nh);

 float **Bayes_matrix(int nrl,int nrh,int ncl,int nch);
 void free_Bayes_matrix(float **m,int nrl,int nrh,int ncl,int nch);
 double **Bayes_dmatrix(int nrl,int nrh,int ncl,int nch);
 void free_Bayes_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);
 int **Bayes_imatrix(int nrl,int nrh,int ncl,int nch);
 void free_Bayes_imatrix(int **m,int nrl,int nrh,int ncl,int nch);

#else

    #define Bayes_vector        vector
    #define free_Bayes_vector   free_vector
	#define Bayes_dvector		dvector
	#define free_Bayes_dvector  free_dvector
    #define Bayes_ivector       ivector
    #define free_Bayes_ivector  free_ivector
    #define Bayes_matrix        matrix
    #define free_Bayes_matrix   free_matrix
	#define Bayes_dmatrix		dmatrix
	#define free_Bayes_dmatrix  free_dmatrix
    #define Bayes_ivector       ivector
    #define free_Bayes_ivector  free_ivector

#endif /* USE_COOLEN_VECTOR */


///////////doesn't really belong here, but required for compilation
//void matrices_error(char error_text[]);
#define matrices_error(error_text) \
  { printf("Matrices library Run-Time Error: %s", (error_text)); \
    return NULL; }

