#include <math.h>

/* (RNG) random number generator from pdf of David Jones, UCL Bioinformatics Group.
  implementation of a RNG proposed by G. Marsaglia. */

static unsigned int x = 125453189, y = 328434370, z = 26488629, c = 1836521; //int x = 123456789, y = 362436000, z = 521288629, c = 7654321; // Seed Variables

unsigned int JKISS(){

	unsigned long long t;

	x = 314527869 * x + 1234567;
	y ^= y << 5; y ^= y >> 7; y ^= y << 22;
	t = 4294584393ULL * z + c; c = t >> 32; z=t;

	return x + y + z; 
}

double uni_dblflt(){
	/*creates a full 53-bit precision double uniform random variable*/
	double x;
	unsigned int a,b;
	
	a = JKISS() >> 6; /* Upper 26 bits */
	b = JKISS() << 5; /* Upper 27 bits */
	x = (a * 134217728.0 + b) / 9007199254740992.0;

	return x;
}

/* Generate gaussian deviate with mean 0 and stdev 1*/

double gaussrnd(){

	double x,y,r;

	do {
		x = 2.0 * uni_dblflt() - 1.0;
		y = 2.0 * uni_dblflt() - 1.0;
		r = x*x + y*y;
	} while(r == 0.0 || r >= 1.0);

	r = sqrt((-2.0 * log(r)) / r);

	return x * r;
}
