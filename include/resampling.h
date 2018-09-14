/*
 * resampling.h
 *
 *  Created on: 5.5.2017
 *      Author: heine
 */

#ifndef RESAMPLING_H_
#define RESAMPLING_H_

double depth_first_sum(double *x, long len);
double *exclusive_prefix_sum(double *w, long len);
double serial_multinomial_resample(long n, double *w, double *w_recv, long *inds);
double serial_multinomial_resample_sng(long n, double *w, long *inds);
double serial_multinomial_resample_noutk(long n, double *w, double* W, long *inds, long k);
double *combining_exclusive_prefix_sum(double *w, double *w_recv, long n);
void set_resamplingseed(int id);
void binary_resampler(long N, double w, long *counts);
double serial_multinomial_resample_sng_noutk(long n, double *w, long *inds, long k) ;
void simplified_resample(long n, long *inds, long k);
void serial_multinomial_resample_counts(long n, long k, double *w, long *cnts);

#endif /* RESAMPLING_H_ */
