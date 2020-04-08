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
