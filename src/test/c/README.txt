SLIMCurve Testing project

I. Test Framework

This test program exercises C code with parameters given in an input file and compares the numerical results against base values.

The input file uses the JSON format.  Here is a sample test from such a file:


        {
            "test":"mono-exp",
            "comment":"Test mono-exp transient 3, noise type 5",
            "tolerance":"0.45",
            "inputs":{
               "noise":"5",
               "xinc":"0.0390625",
               "fitstart":"40",
               "fitstop":"210",
               "ndata":"256",
               "counts":{
                  "size":"256",
                  "value":"0, 3, ... 0, 2 "
               }
            },
            "outputs":{
               "a":"9.673768",
               "t":"12.578522",
               "z":"-3.866590"
            }
         },
        
I've added the ellipsis in the "value" tag, there should actually be 256 values.  (Also note that I made all of my JSON values strings, although JSON does support numerical values.  So '"ndata":"256"' would usually be '"ndata":256'.)

The top level JSON tags here are:

"test" - gives the name of the test to be run
"comment" - description of this specific instance
"tolerance" - how close the numerical results must agree to pass the test, given as a percentage
"inputs"/"outputs" - a list of strings, integers, or floats; also 1D and 2D arrays of integers and floats

There is a tests.c source file that is customized to look for the test name string, grab whatever inputs it needs, run the method under test, and compare results to the outputs.

II.  Running the Tests

There are two test files provided, original.json and test.json.  (They are actually identical except the original.json file has empty output sections [i.e. '"outputs":{}']).

If you run the test with the file original.json (in Windows, '"SLIMCurve Testing" original.json > output') it will run with the old NR-dependent versions of the fitting methods, add the results to the JSON tree, and print out the JSON document to standard output.

You can cut & paste from standard output and build test.json.  Running this version exercises the new code under test and prints out test results.

Unfortunately the JSON document is printed without any whitespace (this is common with JSON).  Fortunately there are many web-based JSON validators available, such as http://jsonformatter.curiousconcept.com/#jsonformatter (or just Google "JSON validator").  It's best to copy and paste the JSON text from the output using original.json into the validator, then copy and paste the validated/formatted results into test.json.

So if you want to change one of the input parameters you have to edit original.json and go through this two stage process building a new test.json.  If you are just tweaking the new code and regression testing just running with the existing test.json will suffice.
