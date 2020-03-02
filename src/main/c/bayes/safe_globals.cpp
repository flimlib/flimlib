#include "safe_globals.h"

#ifdef BAYES_MULTI_THREAD
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
    500.0f, // bayesrapidtauhigh
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
