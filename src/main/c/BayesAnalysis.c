#include "BayesAnalysis.h"
#include <stdio.h>
#include <stdlib.h>

int Bayes_fitting_engine(/* Data in... */
						 float xincr,
						 float *trans,
						 int ndata,
						 int fit_start,
						 int fit_end,
						 float laser_period,
						 float instr[],
						 int ninstr,
						 /* Model... */
						 float param[],
						 int paramfree[],
						 int nparam,
						 /* Data out... */
						 float *fitted,
						 float *residuals,
						 float *error,
						 /* Metadata output */
						 float *minuslogprob,
						 int *nphotons,
						 float *chisq)
{
	int modelType;
	switch (nparam)
	{
	case 3:
		modelType = instr == NULL ? FIT_MONOEXP : FIT_IRFANDMONOEXP;
		break;
	case 5:
		modelType = instr == NULL ? FIT_BIEXP : FIT_IRFANDBIEXP;
		break;
	default:
		// assign unknown if not odd integer >= 7
		modelType = nparam > 7 && ((nparam + 1) & 1) == 1 ? FIT_MULTIEXP : FIT_UNKNOWN;
		break;
	}

#ifdef _abc
	printf_s("xinc: %f, laser_period: %f, ndata: %d,\
		fit_start: %d, fit_end: %d, modeltype: %d,\
		minuslogprob: %f, nphotons: %d, chisq: %f\n",
			 xincr, laser_period, ndata, fit_start, fit_end, modelType, *minuslogprob, *nphotons, *chisq);
	printf_s("trans: ");
	for (int i = 0; i < ndata; i++)
	{
		printf_s("%f, ", trans[i]);
	}
	printf_s("\n");
	printf_s("fitted: ");
	for (int i = 0; i < ndata; i++)
	{
		printf_s("%f, ", fitted[i]);
	}
	printf_s("\n");
	printf_s("residuals: ");
	for (int i = 0; i < ndata; i++)
	{
		printf_s("%f, ", residuals[i]);
	}
	printf_s("\n");
	printf_s("param: ");
	for (int i = 0; i < nparam; i++)
	{
		printf_s("%f, ", param[i]);
	}
	printf_s("\n");
	printf_s("paramfree: ");
	for (int i = 0; i < nparam; i++)
	{
		printf_s("%d, ", paramfree[i]);
	}
	printf_s("\n");
	printf_s("error: ");
	for (int i = 0; i < nparam; i++)
	{
		printf_s("%f, ", error[i]);
	}
	printf_s("\n");
#endif
	if (instr != NULL && ninstr != 0)
	{
		// Setup for fit
		BayesInstrRsp_t IRF = {.paramsfixed = {
								   {.cutofffixed = 0, .delayfixed = 0, .weightfixed = 0, .widthfixed = 0},
								   {.cutofffixed = 0, .delayfixed = 0, .weightfixed = 0, .widthfixed = 0},
								   {.cutofffixed = 0, .delayfixed = 0, .weightfixed = 0, .widthfixed = 0},
							   }};
		BayesIrEstConfig_t IRF_Config = bayes_GetIrEstConfig();

		// // type conversion
		// unsigned* uinstr = malloc(ninstr * sizeof(unsigned));
		// if (!uinstr)
		// 	return BAYES_ERR_PARAM_ESTIMATION_FAILURE;
		// for (int i = 0; i < ninstr; i++)
		// 	uinstr[i] = instr[i];

		int ret = bayes_DoBayesInstrRspEstimation(trans, ndata, xincr, fit_start, fit_end,
											   &nphotons, instr, ninstr, xincr, modelType, &IRF, &IRF_Config,
											   laser_period, param, fitted);
		return ret;
	}
	else
		return bayes_fitting_engine(xincr, laser_period, trans, ndata, fit_start, fit_end,
									param, paramfree, nparam, modelType,
									fitted, residuals, error, minuslogprob, nphotons, chisq);
}
