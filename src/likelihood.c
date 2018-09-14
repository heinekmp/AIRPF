/*
 * likelihood.c
 *
 *  Created on: 5.5.2017
 *      Author: heine
 */
#include <math.h>
#include "typedefs.h"

double likelihood(double x, double y) {
  return -( y - x ) * ( y - x ) / 2 / .25;
  //  return exp(-y * y * exp(-x) / (double) 2 - x / (double) 2);
}

// ****************************************************************
// Evaluate likelihoods
// ****************************************************************
void 
evaluateLikelihoods( 
		    double * x,
		    struct Process process,
		    int state_dim,
		    long data_length,
		    long N,
		    long iteration,
		    double * w,
		    double my_weight )
{
  
  double init_w = log( my_weight ); 
  
  for ( long i = 0; i < N; i++) {
    w[ i ] = init_w;
    
    for ( int d = 0; d < state_dim; d++ ) 
      w[ i ] += likelihood( x[ state_dim * i + d ], 
			    process.y[ data_length * d + iteration ] ) ;
    
    w[ i ] = exp( w[ i ] ) / ( double ) N;
    
  }
  
}

void
normaliseWeights( 
		 double * w,
		 double total_w,
		 long N )
{

  for ( long i = 0; i < N ; i++ ) 
    w[ i ] = w[ i ] / total_w;
    
}
