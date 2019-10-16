%module FLIMLib

%include "arrays_java.i" // arrays (without length)
%include "typemaps.i" // pointers as array

%{
typedef bool boolean;
#define SWIG_JavaArrayInBoolean SWIG_JavaArrayInBool
%}

// macro for 1d array maps
%define ARRMAP(MAPNAME, IS_ARR, HAS_LEN, jType, JType, IO, NULLACC)
// construct the correct typemap name
// TODO: change to readable macros
#undef MAPNAME
#if (IS_ARR == 1)
	#if (HAS_LEN == 1)
		%define MAPNAME (jType ARR_INPUT_LEN_##IO##_##NULLACC[], int ARRLEN) %enddef
	#else
		%define MAPNAME jType ARR_INPUT_LEN_##IO##_##NULLACC[] %enddef
	#endif
#else
	#if (HAS_LEN == 1)
		%define MAPNAME (jType *PTR_INPUT_LEN_##IO##_##NULLACC, int ARRLEN) %enddef
	#else
		%define MAPNAME jType *PTR_INPUT_##IO##_##NULLACC %enddef
	#endif
#endif
%typemap(jstype) MAPNAME "jType[]"
%typemap(javain) MAPNAME "$javainput"
%typemap(javadirectorin) MAPNAME "$jniinput"
%typemap(jtype) MAPNAME "jType[]"
%typemap(jni) MAPNAME %{j##jType##Array%}
%typemap(in) MAPNAME (j##jType *jarr, bool do_clean) %{
#define MAP_1D_ARR_$1_name
	// local reference to the java array and the length of it
	jsize len_$1_name;
	jType *carr_$1_name;

	// $input: jType[](J)
	// $1: jType *ARR_INPUT, $2: int ARRLEN(C)
	if (NULLACC && !$input) {
		// if accepts null input and turns out to be null, don't delete afterwards
		// and skip input array creation
		do_clean = false;
		$1 = NULL;
		len_$1_name = 0;
	}
	else {
		do_clean = true;
		
		len_$1_name = (jsize) JCALL1(GetArrayLength, jenv, $input);
		if (!SWIG_JavaArrayIn##JType##(jenv, &jarr, (jType **)&carr_$1_name, $input))
			exit(1);
		// $*1_type is the type of $1 elements
		// should get optimized with release option
		if (sizeof(jType) != sizeof($*1_type)) {
			// accomondate for e.g. sizeof(bool) != sizeof(int)
			$1 = new $*1_type[len_$1_name];
			for (int i = 0; i < len_$1_name; i++){
				$1[i] = ($*1_type) carr_$1_name[i];
			}
			delete[] carr_$1_name;
		}
		else
			$1 = ($*1_type *) carr_$1_name;
	}

#if (HAS_LEN == 1)
		$2 = len_$1_name;
#endif
%}
%typemap(argout) MAPNAME %{
	if (do_clean$argnum) {
#if (IO != JNI_ABORT)
		// update java array if needed
		int i;
		jsize len = jenv->GetArrayLength($input);
		for (i=0; i<len; i++)
			jarr$argnum[i] = (j##jType)$1[i];
#endif
		jenv->Release##JType##ArrayElements($input, jarr$argnum, IO);
	}
%}
%typemap(freearg) MAPNAME %{
	if (do_clean$argnum) {
		// release the resources before exiting
		delete[] $1;
	}
#undef MAP_1D_ARR_$1_name
%}
%enddef
