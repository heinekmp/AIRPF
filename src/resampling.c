/*
 * resampling.c
 *
 *  Created on: 5.5.2017
 *      Author: heine
 */
#define _XOPEN_SOURCE 600

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include "ziggurat.h"

const gsl_rng_type *T; // Generator type
const gsl_rng *rgen;   // Generator

void set_resamplingseed( int id ) {

  struct timespec spec;

  clock_gettime( CLOCK_REALTIME, &spec );
  
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rgen = gsl_rng_alloc( T );

  gsl_rng_set( rgen, ( unsigned long int ) ( round( spec.tv_nsec / 1.0e6 ) + 223 * id ) );

}

double depth_first_sum(double *x, long len) {
  
  long n = len;
  double *x_tmp, out;
  x_tmp = malloc(sizeof(double) * len);
  for (long i = 0; i < len; i++)
    x_tmp[i] = x[i];
  
  while (n > 1) {
    for (long i = 0; i < n / 2; i++) {
      x_tmp[i] = x_tmp[2 * i] + x_tmp[2 * i + 1];
    }
    if (n % 2 != 0) {
      x_tmp[n / 2] = x_tmp[n - 1];
    }
    n = n / 2 + n % 2;
  }
  out = x_tmp[0];
  free(x_tmp);
  return out;
}

double *depth_first_cumsum(double *x, long len) {
  
  double *W;
  W = malloc(sizeof(double) * len);

  for (long i = 0; i < len; i++) {
    W[i] = depth_first_sum(x, i);
  }
  return W;
}

double *exclusive_prefix_sum(double *w, long len) {
  double *W;
  W = malloc(sizeof(double) * len);
  
  W[0] = 0;
  for (long i = 1; i < len; i++) {
    W[i] = W[i - 1] + w[i - 1];
  }
  
  return W;
}

double *combining_exclusive_prefix_sum(double *w, double *w_recv, long n) {
  
  double *W;
  W = malloc(sizeof(double) * 2 * n);
  
  W[0] = 0;
  
  for (long i = 1; i <= n; i++) {
    W[i] = W[i - 1] + w[i - 1];
  }
  
  for (long i = 0; i < n-1; i++) {
    W[n + i + 1] = W[n + i] + w_recv[i];
  }
  
  return W;
}


double serial_multinomial_resample(long n, double *w, double *w_recv, long *inds) {
  
  double *W;
  double total_weight, lnmax = 0, u;
  long j = 2 * n - 1;
  
  W = combining_exclusive_prefix_sum(w, w_recv, n);
  
  total_weight = log(W[2 * n - 1] + w_recv[n - 1]);
  
  for (long i = n; i >= 1; i--) {

    u = gsl_rng_uniform_pos( rgen );

    lnmax = lnmax + log(u) / (double) i;

    u = exp(total_weight + lnmax);

    while (u < W[j])
      j--;

    inds[i - 1] = j;
  }

  free(W);

  return exp(total_weight);
}

double serial_multinomial_resample_sng(long n, double *w, long *inds) {
  
  double *W;
  double total_weight, lnmax = 0, u;
  long j = n - 1;
  
  W = exclusive_prefix_sum(w, n);
  
  total_weight = log(W[n - 1] + w[n - 1]);
  
  for (long i = n; i >= 1; i--) {

    u = gsl_rng_uniform_pos( rgen ); 

    lnmax = lnmax + log(u) / (double) i;

    u = exp(total_weight + lnmax);

    while (u < W[j])
      j--;

    inds[i - 1] = j;
  }

  free(W);

  return exp(total_weight);
}

double serial_multinomial_resample_sng_noutk(long n, double *w, long *inds, long k) {
  
  double *W;
  double total_weight, lnmax = 0, u;
  long j = k - 1;
  
  W = exclusive_prefix_sum(w, k);
  
  total_weight = log( W[ k - 1 ] + w[ k - 1 ] );
  
  for (long i = n; i >= 1; i--) {

    u = gsl_rng_uniform_pos( rgen ); 

    lnmax = lnmax + log(u) / (double) i;

    u = exp(total_weight + lnmax);

    while (u < W[j])
      j--;

    inds[i - 1] = j;
  }

  free(W);

  return exp(total_weight);
}

double serial_multinomial_resample_noutk(long n, double *w, double* W, long *inds, long k) {
  
  double total_weight, lnmax = 0, u;
  long j = k - 1;
  
  total_weight = log(W[k - 1] + w[k - 1]);
    
  for (long i = n; i >= 1; i--) {

    u = gsl_rng_uniform_pos( rgen ); 

    lnmax = lnmax + log(u) / (double) i;

    u = exp(total_weight + lnmax);

    while (u < W[j])
      j--;

    inds[i - 1] = j;
  }
  
  return -1;
}

void binary_resampler(long n, double w, long *counts)  {

  counts[0] = 0;
  counts[1] = 0;

  for ( long i = 0; i < n; i++ ) 
    counts[  gsl_rng_uniform_pos( rgen ) < w ? 0 : 1 ]++;

}

void 
simplified_resample(
		    long n, 
		    long *inds,
		    long k) 
{

  for ( long i = 0; i < n; i++ ) 
    inds[ i ] = ( long ) ( gsl_rng_uniform_int( rgen, ( unsigned long int ) k ) ); 
}

void serial_multinomial_resample_counts(long n, long k, double *w, long *cnts) {
  
  double *W;
  double total_weight, lnmax = 0, u;
  long j = k - 1;
  
  for( long i = 0; i < k ; i++ ) 
    cnts[ i ] = 0;

  W = exclusive_prefix_sum(w, k);
  
  total_weight = log(W[k - 1] + w[k - 1]);
  
  for (long i = n; i >= 1; i--) {

    u = gsl_rng_uniform_pos( rgen );

    lnmax = lnmax + log(u) / (double) i;

    u = exp(total_weight + lnmax);

    while (u < W[j])
      j--;

    cnts[j]++;
  }

  free(W);
  
}
