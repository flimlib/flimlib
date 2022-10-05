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
#pragma once

#include "bayes_RapidBayesDecayAnalysis.h"
#include "bayes_Types.h"

// thread-local version of global variables in the bayes library
#if __cplusplus
extern "C" {
#endif

BayesRapidValueStore_t* bayes_GetRapidValueStorePtrSafe();

BayesMonoRapidGridConfig_t* bayes_GetMonoRapidGridConfigPtrSafe();

BayesBiRapidGridConfig_t* bayes_GetBiRapidGridConfigPtrSafe();

#if __cplusplus
}
#endif
