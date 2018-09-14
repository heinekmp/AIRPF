#define _XOPEN_SOURCE 600

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <time.h>

void 
setupInitialSeed(
		 int id, 
		 gsl_rng *rgen) 
{

  struct timespec spec;
  
  clock_gettime( CLOCK_REALTIME, &spec );
  
  gsl_rng_set( rgen, (unsigned long int) ( round( spec.tv_nsec / 1.0e6 ) + 20 * id ) );

}

void createInitialSample( 
			 double *x, 
			 long N, 
			 gsl_rng *rgen, 
			 int state_dim ) 
{
  
  for ( long i = 0; i < N; i++)
    for ( int d = 0; d < state_dim; d++ ) 
      x[ state_dim * i + d ] = gsl_ran_gaussian_ziggurat( rgen, 1.0 );
  
}


