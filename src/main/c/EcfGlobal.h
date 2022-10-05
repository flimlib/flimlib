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



#ifdef __cplusplus
extern "C" {
#endif

int GCI_marquardt_global_exps_instr(float xincr, float **trans,
                    int ndata, int ntrans, int fit_start, int fit_end,
                    float instr[], int ninstr,
                    noise_type noise, float sig[], int ftype,
                    float **param, int paramfree[], int nparam,
                    restrain_type restrain, float chisq_delta,
                    float **fitted, float **residuals,
                    float chisq_trans[], float *chisq_global, int *df,
                    int drop_bad_transients);
int GCI_marquardt_global_generic_instr(float xincr, float **trans,
                     int ndata, int ntrans, int fit_start, int fit_end,
                     float instr[], int ninstr,
                     noise_type noise, float sig[],
                     float **param, int paramfree[], int nparam, int gparam[],
                     restrain_type restrain, float chisq_delta,
                     void (*fitfunc)(float, float [], float *, float [], int),
                     float **fitted, float **residuals,
                     float chisq_trans[], float *chisq_global, int *df);

#ifdef __cplusplus
}
#endif
