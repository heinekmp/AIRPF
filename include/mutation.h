/*
 * mutation.h
 *
 *  Created on: 17.5.2017
 *      Author: heine
 */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifndef MUTATION_H_
#define MUTATION_H_

//void set_mutationseed( int id );
void 
mutation( 
	 long N, 
	 double *x, 
	 int state_dim, 
	 gsl_rng * rgen );

#endif /* MUTATION_H_ */
