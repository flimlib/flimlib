#include "safe_globals.h"

#ifdef BAYES_MULTI_THREAD
thread_local BayesRapidValueStore_t bayes_RapidValueStore;
#else
BayesRapidValueStore_t bayes_RapidValueStore;
#endif // BAYES_MULTI_THREAD

// The thread-safe version of the function
BayesRapidValueStore_t* bayes_GetRapidValueStorePtrSafe()
{
	return (&bayes_RapidValueStore);
}

