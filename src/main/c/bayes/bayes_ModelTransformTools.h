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
#ifndef BAYES_MODEL_TRANSFORM_TOOLS_H
#define BAYES_MODEL_TRANSFORM_TOOLS_H


/* Parameter value conversion between Bayesian and conventional models... */

float bayes_ToBayesModelTransformFromParamZ(float Z, int nbins, int nphotons);

float bayes_FromBayesModelTransformToParamZ(float w0, int nbins, int nphotons);

float bayes_FromBayesModelWeightAndTauToParamA( float w_i,
                                                float tau_i,
                                                float binsize,
                                                float interval,
                                                float delay,
                                                int   nphotons);

float bayes_ToBayesModelWeightFromParamAAndTau( float A_i,
                                                float tau_i,
                                                float binsize,
                                                float interval,
                                                float delay,
                                                int   nphotons);

int bayes_CheckForValidBayesModelWeightParamValues(/* Model and params... */
                                                   int     nweights,
                                                   double *weights,
                                                   int    *isweightfixed);

int bayes_CheckForValidBayesModelLifetimeParamValues(/* Model and params... */
                                                     int     nlifetimes,
                                                     double *lifetimes,
                                                     int    *islifetimefixed);

int bayes_ToBayesModelParamValuesFromConventionalModelParamValues(/* Model and params... */
                                                                  int    nparams,
                                                                  float *params_bayes,
                                                                  float *params_conventional,
                                                                  /* Instrument response... */
                                                                  float  delay,
                                                                  /* Transient... */
                                                                  int    nbins,
                                                                  int    nphotons,
                                                                  float  interval);

int bayes_ToConventionalModelParamValuesFromBayesModelParamValues(/* Model and params... */
                                                                  int    nparams,
                                                                  float *params_bayes,
                                                                  float *params_conventional,
                                                                  /* Instrument response... */
                                                                  float  delay,
                                                                  /* Transient... */
                                                                  int    nbins,
                                                                  int    nphotons,
                                                                  float  interval);

int bayes_ToConventionalModelParamValuesFromBayesModelParamValuesFreeParamsOnly(/* Model and params... */
                                                                                int    nparams,
                                                                                int   *isweightfixed,
                                                                                int   *istaufixed,
                                                                                float *params_bayes,
                                                                                float *params_conventional,
                                                                                /* Instrument response... */
                                                                                float  delay,
                                                                                /* Transient... */
                                                                                int    nbins,
                                                                                int    nphotons,
                                                                                float  interval);

/* Weight and lifetime extraction utilities... */

int bayes_AllocateWeightsAndTausVectors(int nparams, int *ndecays, double **weights, double **taus);

int bayes_AllocateWeightsAndTausIsFixedVectors(int ndecays, int **isweightfixed, int **istaufixed);

int bayes_PopulateWeightsAndTausVectorsFromParamVector(int nparams, float *params, double *weights, double *taus);

int bayes_PopulateParamVectorFromWeightsAndTausVectors(float *params, int nparams, double *weights, double *taus);

void bayes_FreeWeightsAndTausVectors(int ndecays, double *weights, double *taus);

void bayes_FreeWeightsAndTausFixedVectors(int ndecays, int *isweightfree, int *istaufree);


/* Free and fixed parameter utilities... */

#define BAYES_PARAM_VALUE_FREE           0
#define BAYES_PARAM_VALUE_USER_FIXED     1
#define BAYES_PARAM_VALUE_MODEL_FIXED    2

int bayes_DetermineModelParamsFixed(int nparams, int *isparamuserfixed, int *nparamsfixed, int *isparamfixed, int ndecays,                                    
                                    int *isweightfixed, int *nweightsfixed, int *istaufixed, int *ntausfixed);

int bayes_AllocateFreeAndFixedParamVectors(int nxfree, double **xfree, int nxfixed, double **xfixed);

int bayes_PopulateParamVectorFromFreeAndFixedVectors(float *x, int nx, double *xfree, int nxfree,
                                                     double *xfixed, int nxfixed, int *isfixed);

int bayes_PopulateFreeAndFixedVectorsFromParamVector(float *x, int nx, double *xfree, int nxfree,
                                                     double *xfixed, int nxfixed, int *isfixed);

int bayes_UpdateWeightsVectorModelDefinedValue(double *weights, int nweights, int *isweightfixed);

void bayes_FreeFreeAndFixedParamVectors(int nxfree, double *xfree, int nxfixed, double *xfixed);



/* Parameter random initialization... */

int bayes_RandomlyInitWeightsVector(double *weights, int nweights, int nweightsfixed, int *isweightfixed);

int bayes_RandomlyInitTausVector(double *taus, int ntaus, int *istaufixed, double interval);



/* Miscellaneous model parameter related tools... */

int bayes_OrderDecaysByDecreasingLifetimes(int ndecays, double *weights, double *taus);




float bayes_ComputeRawBackgroundZ(float w0,
                                  float binsize,
                                  float interval,
                                  int nphotons);


float bayes_ComputeRawAmplitudeA(float w0,
                                 float w1,
                                 float binsize,
                                 float interval,
                                 float delay,
                                 int nphotons);

#if 0
float bayes_ModelTransformComputePhotonCount(float delay,
			                                 float width,
				                             float interval,
				                             float t,
				                             float w0,
				                             float w1,
				                             float binsize,
				                             int   nphotons);
#else
float bayes_ModelTransformComputePhotonCount(float interval,
				                             float w0,
                                             float likelihood,
				                             float binsize,
				                             int   nphotons);
#endif

//void scaleDataAccordingToSignalNoise (float *data, int n, float *signal);


#endif /* BAYES_MODEL_TRANSFORM_TOOLS_H */
