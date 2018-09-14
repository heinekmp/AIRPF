/*
 * likelihood.h
 *
 *  Created on: 5.5.2017
 *      Author: heine
 */
#include "typedefs.h"

#ifndef LIKELIHOOD_H_
#define LIKELIHOOD_H_

double 
likelihood(
	   double x, 
	   double y);

void 
evaluateLikelihoods( 
		    double * x,
		    struct Process process,
		    int state_dim,
		    long data_length,
		    long N,
		    long iteration,
		    double * w,
		    double my_weight );
void
normaliseWeights( 
		 double * w,
		 double total_w,
		 long N );


#endif /* LIKELIHOOD_H_ */
