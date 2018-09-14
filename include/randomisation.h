#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void 
setupInitialSeed( 
		 int id, 
		 gsl_rng *rgen );

void createInitialSample( 
			 double *x, 
			 long N, 
			 gsl_rng *rgen, 
			 int state_dim );
