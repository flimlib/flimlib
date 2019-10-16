#include "EcfInternal.h"
#include "json.h"
#include "parser.h"
#include "output.h"
#include "tests.h"
#include <assert.h>
#include <locale.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#define TRIALS 10000000

//TODO ARG this is written in a crap style that Visual Studio 2008 will not like; i.e. defining variables anywhere
//TODO ARG after mono-exp and tri-exp are working, port to VS2008

char *report_noise_type(noise_type noise) {
	char *return_value;
	switch (noise) {
		case NOISE_CONST:
			return_value = "NOISE_CONST";
			break;
		case NOISE_GIVEN:
			return_value = "NOISE_GIVEN";
			break;
		case NOISE_POISSON_DATA:
			return_value = "NOISE_POISSON_DATA";
			break;
		case NOISE_POISSON_FIT:
			return_value = "NOISE_POISSON_FIT";
			break;
		case NOISE_GAUSSIAN_FIT:
			return_value = "NOISE_GAUSSIAN_FIT";
			break;
		case NOISE_MLE:
			return_value = "NOISE_MLE";
			break;
	}
	return return_value;
}

char *report_fit_type(fit_type fit) {
	char *return_value;
	switch (fit) {
		case FIT_GLOBAL_MULTIEXP:
			return_value = "FIT_GLOBAL_MULTIEXP";
			break;
		case FIT_GLOBAL_STRETCHEDEXP:
			return_value = "FIT_GLOBAL_STRETCHEDEXP";
			break;
	}
	return return_value;
}

char *report_restrain_type(restrain_type restrain) {
	char *return_value;
	switch (restrain) {
		case ECF_RESTRAIN_DEFAULT:
			return_value = "ECF_RESTRAIN_DEFAULT";
			break;
		case ECF_RESTRAIN_USER:
			return_value = "ECF_RESTRAIN_USER";
			break;
	}
	return return_value;
}

double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks*10)/CLOCKS_PER_SEC;
	return diffms;
}

float **make_copy_2D(float data[3][3], int n) {
	float **copy = GCI_ecf_matrix(n, n);
	int i, j;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			copy[i][j] = data[i][j];
		}
	}
	return copy;
}

float *make_copy_1D(float data[3], int n) {
	int size = n * sizeof(float);
	float *copy = (float *) malloc(size);
	int i;
	for (i = 0; i < n; ++i) {
		copy[i] = data[i];
	}
	return copy;
}

void compare_solve_timings() {
	clock_t begin, end;
	int trial;
	float a[3][3] = {
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0 };
	float b[3] = {
		1.0, 2.0, 3.0 };
	
    // insert code here...
    printf("FLIMLib testing\n");
		
	begin = clock();
	for (trial = 0; trial < TRIALS; ++trial) {
        float **tmp_a= make_copy_2D(a, 3);
		float *tmp_b = make_copy_1D(b, 3);
		
		GCI_solve_Gaussian(tmp_a, 3, b);
		
		GCI_ecf_free_matrix(tmp_a);
		free(tmp_b);
	}
	end = clock();
	printf("time %f\n", diffclock(end, begin));
	
	
	begin = clock();
	for (trial = 0; trial < TRIALS; ++trial) {
        float **tmp_a= make_copy_2D(a, 3);
		float *tmp_b = make_copy_1D(b, 3);
		
		GCI_solve_lu_decomp(tmp_a, 3, b);
		
		GCI_ecf_free_matrix(tmp_a);
		free(tmp_b);    }
	end = clock();
	printf("time %f\n", diffclock(end, begin));
	
    return;
}

/*
 * Processes a single test.
 */
void process(json_t *cursor, float tolerance, int *changed, int *success)
{
	char *test;
	char *comment;
	json_t *inputs;
	json_t *outputs;

	test = getString(cursor, "test");	
	inputs = getNamedChildValue(cursor, "inputs");
	outputs = getNamedChildValue(cursor, "outputs");
	
	printf("\nTEST: %s\n", test);
	comment = getString(cursor, "comment");
	if (NULL != comment) {
		printf("%s\n", comment);
	}
	tolerance = getFloat(cursor, "tolerance", tolerance);
	
	if (NULL != test && NULL != inputs)
	{
		do_test(test, inputs, outputs, tolerance, changed, success);
	}
}

/*
 * Main program.  Reads in a JSON test file and runs the tests specified within.
 */
int main (int argc, char **argv)
{
	int changed = FALSE;
	int success = TRUE;
	int tests = 0;
	int successes = 0;
	json_t *document = NULL;
	json_t *cursor = NULL;
	FILE *fp;
	int i;
	float tolerance;
	
	if (argc < 2)
	{
		printf("usage: test document1.json ...\n");
		return EXIT_SUCCESS;
	}
	
	for (i = 1; i < argc; i++)
	{
		printf("PROCESSING FILE %s...\n",argv[i]);
		fp = fopen(argv[i],"r");
		if(fp == NULL)
		{
			printf("file \"%s\" couldn't be opened\n",argv[i]);
		}
		else
		{
			switch(json_stream_parse(fp, &document))
			{
				case JSON_OK:
					cursor = document;
					assert(cursor->type == JSON_OBJECT);
					cursor = cursor->child;
					assert(cursor->type == JSON_STRING);
					if (strcmp("testrun", cursor->text) != 0) {
						printf("Missing 'testrun' element\n");
					}
					cursor = cursor->child;
					tolerance = getFloat(cursor, "tolerance", 0.0f);
					if (0.0f == tolerance) {
						printf("Missing 'tolerance' element of 'testrun', using 10 percent\n");
						tolerance = 10.0f; // within 10% agreement of results
					}
					
					// loop over tests
					cursor = getNamedChildValue(cursor, "tests");
					if (NULL != cursor && JSON_ARRAY == cursor->type) {
						cursor = cursor->child;
						while (NULL != cursor) {
							success = TRUE;
							process(cursor, tolerance, &changed, &success);
							cursor = cursor->next;
							
							// account for this test
							++tests;
							if (success) {
								++successes;
							}
						}
					}

					// if the JSON document has changed (test outputs added) write it out
					if (changed) {
						printf("JSON document changed:\n\n");
						json_stream_output(stdout,document);
						printf("\n\n");
					}
					else {
						printf("\n%d TESTS RUN, %d SUCCESSES, %d FAILURES\n\n", tests, successes, tests - successes);
					}
					
					json_free_value(&document);
					break;
					
				default:
					printf("Problems parsing JSON input\n");
					break;
			}
		}
		fclose(fp);
	}
	
	//compare_solve_timings();

	return EXIT_SUCCESS;
}


