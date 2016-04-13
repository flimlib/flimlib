/* File: EcfGlobal.i */
%module EcfGlobal 
/* For things like float* */
/*%include "typemaps.i"*/
%include "cpointer.i"
%include "arrays_java.i"

%{
#include "Ecf.h"
#include "EcfGlobal.h"
%}

%pointer_functions(float, floatp);
%pointer_functions(int, intp);
%javaconst(1);
%include "Ecf.h"
%include "EcfGlobal.h"
