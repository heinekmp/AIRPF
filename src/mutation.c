/*
 * mutation.c
 *
 *  Created on: 5.5.2017
 *      Author: heine
 */

#define _XOPEN_SOURCE 600
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <mpi.h>

//static const gsl_rng_type * T; // Generator type
//static const gsl_rng * rgen;   // Generator

/*void set_mutationseed( int id ) {

  struct timespec spec;

  clock_gettime( CLOCK_REALTIME, &spec );
  
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rgen = gsl_rng_alloc( T );

  gsl_rng_set( rgen, ( unsigned long int ) ( round( spec.tv_nsec / 1.0e6 ) + 34 * id ) );
  }*/

void 
mutation( 
	 long N, 
	 double *x, 
	 int state_dim,
	 gsl_rng * rgen ) 
{
  
  for ( long i = 0; i < N; i++ ) 
    for ( int d = 0; d  < state_dim; d++ ) 
      x[ state_dim *  i + d ] = x[ state_dim * i + d ] + 
	gsl_ran_gaussian_ziggurat( rgen, 1.0 );

}
