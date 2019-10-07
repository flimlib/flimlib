#pragma once

#include "bayes_RapidBayesDecayAnalysis.h"

// thread-local version of global variables in the bayes library
#if __cplusplus
extern "C" {
#endif

BayesRapidValueStore_t* bayes_GetRapidValueStorePtrSafe();

#if __cplusplus
}
#endif
