%module FLIMLib

%include "enums.swg"

%define ENUMMAP(MAPNAME, ctype, JType)
%rename(JType) ctype;
#undef MAPNAME
// Conversion: JType(J) -> int(C) in arguments
%define MAPNAME int ENUM_NAME %enddef
%typemap(jstype) MAPNAME "JType"
%typemap(javain) MAPNAME "$javainput.swigValue()"
%typemap(jtype) MAPNAME "int"
%typemap(jni) MAPNAME "jint"
%typemap(in) MAPNAME {
	$1 = $input;
}

%enddef
