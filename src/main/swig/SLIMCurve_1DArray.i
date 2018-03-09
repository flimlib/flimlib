%module SLIMCurve

%include "arrays_java.i" // arrays (without length)
%include "typemaps.i" // pointers as array

// macro for 1d array maps, IO specifies how java array is updated
%define ARRMAP(mapName, jType, JType, IO)
	// Conversion: jType[](J) -> (jType [], int)(C) for input arguments
%define mapName (jType ARR_INPUT[], int ARRLEN) %enddef
%typemap(jstype) mapName "jType[]"
%typemap(javain) mapName "$javainput"
%typemap(javaout) mapName {return $jnicall;}
%typemap(jtype) mapName "jType[]"
	// cannot use macro inside " "
%typemap(jni) mapName "JXXXARRAY(jType)"
%typemap(in) mapName (j##jType *jarr) {
	
	// $input: jType[](J)
	// $1: jType *ARR_INPUT, $2: int ARRLEN(C)
	if (!SWIG_JavaArrayIn##JType##(jenv, &jarr, (jType **)&$1, $input)) return 0;
	$2 = (int) JCALL1(GetArrayLength, jenv, $input);
}
%typemap(freearg) mapName {
	// release the resources before exiting
	delete[] $1;
}
%enddef
