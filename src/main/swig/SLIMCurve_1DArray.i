%module SLIMCurve

%include "arrays_java.i" // arrays (without length)
%include "typemaps.i" // pointers as array

// macro for 1d array maps
%define ARRMAP(MAPNAME, IS_ARR, HAS_LEN, jType, JType, IO, NULLACC, KEEP_ALLOC)
// construct the correct typemap name
// TODO: change to readable macros
#undef MAPNAME
#if IS_ARR == 1
	#if HAS_LEN == 1
		%define MAPNAME (jType ARR_INPUT[], int ARRLEN) %enddef
	#else
		%define MAPNAME jType ARR_INPUT[] %enddef
	#endif
#else
	#if HAS_LEN == 1
		%define MAPNAME (jType *ARR_INPUT, int ARRLEN) %enddef
	#else
		%define MAPNAME jType *ARR_INPUT %enddef
	#endif
#endif
%typemap(jstype) MAPNAME "jType[]"
%typemap(javain) MAPNAME "$javainput"
%typemap(javadirectorin) MAPNAME "$jniinput"
%typemap(jtype) MAPNAME "jType[]"
	// cannot use ## macro inside " "
%typemap(jni) MAPNAME "JXXXARRAY(jType)"
%typemap(in) MAPNAME (j##jType *jarr, void *combined_arr = 0, bool do_clean) {
	// local reference to the java array and the length of it
	jsize len;
	jType *carr;

	// $input: jType[](J)
	// $1: jType *ARR_INPUT, $2: int ARRLEN(C)
	if (NULLACC && !$input) {
		// if accepts null input and turns out to be null, don't delete afterwards
		// and skip input array creation
		do_clean = false;
		$1 = NULL;
		len = 0;
	}
	else {
		do_clean = true;
		
		len = (jsize) JCALL1(GetArrayLength, jenv, $input);
		if (!SWIG_JavaArrayIn##JType##(jenv, &jarr, (jType **)&carr, $input))
			exit(1);
#if KEEP_ALLOC == 1
// "KEEP_ALLOC" is used to keep the array alive in case it is used later
// (e.g. allocated for a DecayModelSelParamValuesAndFit struct)
		// make sure there is enough space for the jarray pointer
		combined_arr = std::malloc(len * sizeof(jType) + sizeof(void*));
		combined_arr = &$input;
		// skip the pointer
		$1 = (jType*)((char*)combined_arr + sizeof(void*));
		for (int i = 0; i < len; i++)
			$1[i] = (j##jType) carr[i];
		printf("%p %p\n", combined_arr, $1);
#else
		$1 = carr;
#endif
	}
#if HAS_LEN == 1
		$2 = len;
#endif
}
%typemap(argout) MAPNAME %{
	if (do_clean$argnum) {
		// update java array if the array is not needed any more
		int i;
		jsize len = jenv->GetArrayLength($input);
		for (i=0; i<len; i++)
			jarr$argnum[i] = (j##jType)$1[i];
		jenv->Release##JType##ArrayElements($input, jarr$argnum, IO);
	}
%}
%typemap(freearg) MAPNAME %{
	if (do_clean$argnum) {
		// release the resources before exiting
#if KEEP_ALLOC == 1
		std::free(combined_arr$argnum);
#else
		delete[] $1;
#endif
	}
%}
%enddef
