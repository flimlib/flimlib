slim-curve-test
===============

SLIM Curve test program

This is a program that exercises SLIM Curve fitting routines based on data from a json input file.  The program runs in two modes:  Usually the json file has a series of tests with expected results and tolerances for each test against which the calculated results are compared, so that the test either passes or fails.  In the other mode the json file just specifies a series of tests and a new json file gets created with the calculated results saved.  This second mode is a means of creating a set of baseline results.

The slim-curve main repository is required alongside this slim-curve-test repository.
