#include "BayesAnalysis.h"
#include <stdio.h>
#include <stdlib.h>

static const float RANGE_FACTOR = 2;

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
		BayesMonoRapidGridConfig_t *defaultMonoGridConfig = bayes_GetMonoRapidGridConfigPtrSafe();

		defaultMonoGridConfig->bayesrapidwlow = param[1] / RANGE_FACTOR;
		defaultMonoGridConfig->bayesrapidtaulow = param[2] / RANGE_FACTOR;

		defaultMonoGridConfig->bayesrapidbghigh = param[0] * RANGE_FACTOR;
		defaultMonoGridConfig->bayesrapidwhigh = param[1] * RANGE_FACTOR;
		defaultMonoGridConfig->bayesrapidtauhigh = param[2] * RANGE_FACTOR;
		break;
	case 5:
		modelType = instr == NULL ? FIT_BIEXP : FIT_IRFANDBIEXP;
		BayesBiRapidGridConfig_t *defaultBiGridConfig = bayes_GetBiRapidGridConfigPtrSafe();
		defaultBiGridConfig->bayesrapidbibgmin = param[0] / RANGE_FACTOR;
		defaultBiGridConfig->bayesrapidbiw1low = param[1] / RANGE_FACTOR;
		defaultBiGridConfig->bayesrapidbitau1low = param[2] / RANGE_FACTOR;
		defaultBiGridConfig->bayesrapidbiw2low = param[3] / RANGE_FACTOR;
		defaultBiGridConfig->bayesrapidbitau2low = param[4] / RANGE_FACTOR;

		defaultBiGridConfig->bayesrapidbibgmax = param[0] * RANGE_FACTOR;
		defaultBiGridConfig->bayesrapidbiw1high = param[1] * RANGE_FACTOR;
		defaultBiGridConfig->bayesrapidbitau1high = param[2] * RANGE_FACTOR;
		defaultBiGridConfig->bayesrapidbiw2high = param[3] * RANGE_FACTOR;
		defaultBiGridConfig->bayesrapidbitau2high = param[4] * RANGE_FACTOR;
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

		int ret = bayes_DoBayesInstrRspEstimation(trans, ndata, xincr, fit_start, fit_end,
												  nphotons, instr, ninstr, xincr, modelType, &IRF, &IRF_Config,
												  laser_period, param, fitted);
		return ret;
	}
	else
		return bayes_fitting_engine(xincr, laser_period, trans, ndata, fit_start, fit_end,
									param, paramfree, nparam, modelType,
									fitted, residuals, error, minuslogprob, nphotons, chisq);
}
