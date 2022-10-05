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
#include "safe_globals.h"

#ifndef BAYES_SINGLE_THREAD
#define SAFE_GLOBAL thread_local
#else
#define SAFE_GLOBAL
#endif

static SAFE_GLOBAL BayesRapidValueStore_t bayes_RapidValueStore;

// Search constraints
static SAFE_GLOBAL BayesMonoRapidGridConfig_t bayes_MonoRapidGridConfig = {
    // c++(<20) doesn't accept designated initializer list
    100,    // bayesrapidtaupts
    1e-6f,  // bayesrapidtaulow
    5.0f, // bayesrapidtauhigh
    200,    // bayesrapidwpts
    0.0f,   // bayesrapidwlow
    1.0f,   // bayesrapidwhigh
    10.0f   // bayesrapidbghigh
};

static SAFE_GLOBAL BayesBiRapidGridConfig_t bayes_BiRapidGridConfig = {
    100,       // bayesrapidbitaupts
    0.01f,     // bayesrapidbitaulow
    4.0f,      // bayesrapidbitauhigh
    50,        // bayesrapidbiweightpts
    0.0f,      // bayesrapidbiweightlow
    1.0f,      // bayesrapidbiweighthigh
    0.0f,      // bayesrapidbibgmin
    0.244898f, // bayesrapidbibgmax

    0.0f,      // bayesrapidbiw0low
    0.489796f, // bayesrapidbiw1low
    0.0f,      // bayesrapidbiw2low
    1.501212f, // bayesrapidbitau1low
    0.050303f, // bayesrapidbitau2low

    0.102041f, // bayesrapidbiw0high
    1.0f,      // bayesrapidbiw1high
    0.489796f, // bayesrapidbiw2high
    2.508788f, // bayesrapidbitau1high
    1.017576f  // bayesrapidbitau2high
};

// The thread-safe version of the functions
BayesRapidValueStore_t *bayes_GetRapidValueStorePtrSafe()
{
    return (&bayes_RapidValueStore);
}

BayesMonoRapidGridConfig_t *bayes_GetMonoRapidGridConfigPtrSafe()
{
    return (&bayes_MonoRapidGridConfig);
}

BayesBiRapidGridConfig_t *bayes_GetBiRapidGridConfigPtrSafe()
{
    return (&bayes_BiRapidGridConfig);
}
