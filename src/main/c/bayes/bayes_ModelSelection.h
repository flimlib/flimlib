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
