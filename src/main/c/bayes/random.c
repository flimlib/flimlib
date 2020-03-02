
// CVI Random function seems to be better than the rand1 and rand2 functions that were used
// according to the dieharder tests: http://www.phy.duke.edu/~rgb/General/dieharder.php
// See labnotes for 16th and 17th March 2017
// I now use my implementation of the Wichmann-Hill algorithm
// P Barber, March 2017
// Wichmann, Brian; Hill, David (1982). "Algorithm AS 183: An Efficient and Portable Pseudo-Random Number Generator". Journal of the Royal Statistical Society. Series C (Applied Statistics), Vol. 31, No. 2 (1982), pp. 188-190.

#include <stdlib.h>
#include <math.h>
#include <time.h>

static unsigned int wh_IX, wh_IY, wh_IZ;  

static void seed_WichmannHill(unsigned long seed)
{
	unsigned int randomSeed;           
    
	if (seed)
        randomSeed = seed; 
    else
	    do
	    {
	        randomSeed = clock();
	    } while (!randomSeed);    


	srand(randomSeed);
	wh_IX = rand();
	wh_IY = rand();
	wh_IZ = rand();

}

static double rand_WichmannHill(void)
{
    double randNum;

    wh_IX = (171 * wh_IX) % 30269;
    wh_IY = (172 * wh_IY) % 30307;
    wh_IZ = (170 * wh_IZ) % 30323;

    randNum = wh_IX / 30269.0 + wh_IY / 30307.0 + wh_IZ / 30323.0;
    randNum -= floor(randNum);
    if(randNum < 0) randNum += 1.0; // make sure 0-1
	
    return randNum;
}

void rand_InitializeRandomSeed(void)
{
	// pass zero to set by clock
	seed_WichmannHill(0);
}

float rand_RandomFloat(void)
{
	return((float)rand_WichmannHill());
}

double rand_RandomDouble(void)
{
	return(rand_WichmannHill());
}



