#ifndef BAYES_MODEL_SELECTION_H
#define BAYES_MODEL_SELECTION_H


#include "bayes_Types.h"
#include "bayes_Interface.h"


int bayes_DetemineDecayModelRelativeLikelihoods(
	/* Data in... */
	int                       nbins,
	int                       fitstart,
	double                   *binwalls,
	int                      *data,
	int                       nphotons,
	/* Instrumentation... */
	double                    interval,
	double                    modperiod,
	BayesInstrRsp_t          *instr,
	/* RLD-derived estimates for search initialisation... */
	BayesParamValsAndFit_t   *decayestimates,
	/* Model selection data out... */
	float                    *decaymodellikelihoods,
	BayesParamValsAndFit_t   *paramvalsandfits,
	/* Configuration... */
	int                       rapidanalysis,
	BayesRapidValueStore_t   *rapidgrid);

int bayes_DetemineDecayModelEvidence(
	int              ndecays,
	double          *weights,
	double          *taus,
	double          *hyperparams,
	double           minuslogprob,
	int              nbins,
	double          *binwalls,
	int             *data,
	double           interval,
	double           modperiod,
	BayesInstrRsp_t *instr,
	double          *logmodelevidence);

#endif /* BAYES_MODEL_SELECTION_H */