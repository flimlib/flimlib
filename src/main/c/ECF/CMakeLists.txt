CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(SLIM_CURVE)

SET(ECF_SRC		 	EcfGlobal.c
 					EcfSingle.c
 					EcfSPA.c
 					EcfUtil.c
					ecf.h 
					EcfInternal.h
)

SOURCE_GROUP(slim-curve FILES ${ECF_SRC})

ADD_LIBRARY(slim-curve ${ECF_SRC})

#ADD_EXECUTABLE(slim-curve_test FILES test.c ${ECF_SRC})

