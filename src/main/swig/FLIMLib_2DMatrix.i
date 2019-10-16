%module FLIMLib

%{
#include "ParamMatrix.h"

#define JXXXARRAY(xxx) j##xxx##Array
#define XXX2DMATRIX(xxx) xxx##2DMatrix
typedef ParamMatrix<float> Float2DMatrix;
typedef ParamMatrix<int> Int2DMatrix;
%}

%extend ParamMatrix {
	ParamMatrix *asArray() {
		return self;
	}
}

/* Xxx2DMatrix */
%define MATMAP(mapName, jType, JType, JSignature, MatName)
// Conversion: XXX2DMATRIX(JType)(J) -> jType**(C) in arguments
%typemap(jstype) jType** "MatName"
%typemap(javain) jType** "$javainput.getCPtr($javainput)"
%typemap(jtype) jType** "long"
%typemap(jni) jType** "jlong"
%typemap(in) jType** {
	$1 = ((XXX2DMATRIX(JType)*)$input)->arr;
}
// Conversion: jType[][](J) -> Xxx2DMatrix(C)
%typemap(jstype) (jType **arr, int row, int col) "jType[][]"
%typemap(javain, pre="checkArray($javainput);") (jType **arr, int row, int col) "$javainput"
%typemap(jtype) (jType **arr, int row, int col) "jType[][]"
%typemap(jni) (jType **arr, int row, int col) "jobjectArray"
	// for this special set of arguments (XXX2DMATRIX(JType)::XXX2DMATRIX(JType))
%typemap(in) (jType **arr, int row, int col) {
	// $input: jType[][](J)
	// $1: jType **arr, $2: int row, $3: int col(C)
	$2 = JCALL1(GetArrayLength, jenv, $input);
	JXXXARRAY(jType) arow = (JXXXARRAY(jType))JCALL2(GetObjectArrayElement, jenv, $input, 0);
	// assume rectangular array, else not matrix
	$3 = JCALL1(GetArrayLength, jenv, arow);
	$1 = ParamMatrix<jType>::new_mat($2, $3);
	for (int i = 0; i < $2; i++) {
		JXXXARRAY(jType) arow = (JXXXARRAY(jType))JCALL2(GetObjectArrayElement, jenv, $input, i);
		j##jType* elements = JCALL2(Get##JType##ArrayElements, jenv, arow, 0);
		for (int j = 0; j < $3; j++) {
			$1[i][j] = elements[j];
		}
		// release elements (alloc'ed in Get##JType##ArrayElement)
		// without copying changes back to jType $2 in JVM memory
		JCALL3(Release##JType##ArrayElements, jenv, arow, elements, JNI_ABORT);
	}
}
// Conversion: jType[][](J) <- Xxx2DMatrix(C) in Xxx2DMatrix::asArray only
%typemap(jstype) ParamMatrix<jType> *asArray "jType[][]"
%typemap(javaout) ParamMatrix<jType> *asArray {return $jnicall;}
%typemap(jtype) ParamMatrix<jType> *asArray "jType[][]"
%typemap(jni) ParamMatrix<jType> *asArray "jobjectArray"
%typemap(out) ParamMatrix<jType> *asArray {
	// $1: Xxx2DMatrix(C)
	// $result: jType[][](J)
	// This map just wraps the internal ** around a jType[][]
	int row = $1->nrow;
	int col = $1->ncol;
	const j##jType **data = (const j##jType **)$1->arr;
	jclass _2MClass = JCALL1(FindClass, jenv, "[JSignature");
	$result = JCALL3(NewObjectArray, jenv, row, _2MClass, NULL);
	for (int i = 0; i < row; i++) {
		JXXXARRAY(jType) arow = JCALL1(New##JType##Array, jenv, col);
		JCALL4(Set##JType##ArrayRegion, jenv, arow, 0, col, data[i]);
		JCALL3(SetObjectArrayElement, jenv, $result, i, arow);
	}
}
// Conversion: XXX2DMATRIX(JType)(J) -> (jType **, int, int)(C) for output arguments
%define mapName (jType **ARR_INPUT, int NCOL, int NROW) %enddef
%typemap(jstype) mapName "MatName"
%typemap(javain) mapName "$javainput.getCPtr($javainput)"
%typemap(jtype) mapName "long"
%typemap(jni) mapName "jlong"
%typemap(in) mapName {
	// $input: jType[](J)
	// $1: jType **, $2: int NCOL, $3: int NROW(C)
	XXX2DMATRIX(JType) *fp = (XXX2DMATRIX(JType)*)$input;
	$1 = fp->arr;
	$2 = fp->ncol;
	$3 = fp->nrow;
}
/* Java code to be inserted into FLIMLib class */
%typemap(javaimports) ParamMatrix<jType> %{
import java.util.Arrays;
%}
%typemap(javacode) ParamMatrix<jType> %{
/* Checks array dimensions before converted into a Float2DMatrix */
private static void checkArray(final jType[][] arr) {
	if (arr == null)
		throw new NullPointerException("Array is null");
	final int row = arr.length;
	if (row == 0)
		throw new IllegalArgumentException("Array should have at least 1 row");
	final int col = arr[0].length;
	if (col == 0)
		throw new IllegalArgumentException("Array should have at least 1 column");
	for (jType[] arrRow : arr)
		if (arrRow.length != col)
			throw new IllegalArgumentException("Array not rectangular");
}

@Override
public String toString() {
	// nparray-like string reprensentation
	String data = "";
	jType[][] arr = this.asArray();
	for(int i = 0; i < arr.length; i++) {
		String row = "";
		for(int j = 0; j < arr[i].length; j++) {
			jType abs = Math.abs(arr[i][j]);
			row += String.format((abs >= 1e6 || abs <= 1e-6) ? "%13e " : "%13f ", arr[i][j]);
		}
		if (i == 0)
			data += "[ [ " + row + "]\n";
		else if (i < arr.length - 1)
			data += "  [ " + row + "]\n";
		else
			data += "  [ " + row + "] ]\n";
	}
	return data;
}
%}
%template(JType##2DMatrix) ParamMatrix<jType>;
%enddef

// Use XX2DMatrix instead
%ignore GCI_ecf_matrix;
%ignore GCI_ecf_free_matrix;
%ignore ParamMatrix::arr;
%ignore ParamMatrix::new_mat;
%include "../cpp/ParamMatrix.h"
