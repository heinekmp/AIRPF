/*
 * fileio.h
 *
 *  Created on: 17.5.2017
 *      Author: heine
 */
#include "typedefs.h"

#ifndef FILEIO_H_
#define FILEIO_H_

void read_true_normalisers(double *mu);
void read_true_normaliser_vars(double *sd);
void reset_output();
void output(double *w, double *x, long N, int id, long k, int s);
void writeIndices(long *inds, double *w, long N);
//void writeProcessFile(struct Process process, int id);
void waitFilter(int id, long k, int s);
void combineSamples(double *like, double *x, long N, int id, long k, int stage, double *comb_like,
		double *comb_x);
long countlines(char *filename);
void removeOldFiles(long k, int nid, int S);
struct Estimates estimateOutput( double *w, double *x, int id, long N, int rund_idx, int state_dim );
struct Process readProcessFile( int state_dim );
void combinedEstimateOutput(double x, double w, int run_idx, int world_size, long N);
void writeOutputFile(double *data, long data_length, int world_size, int run_idx, long N, short intact, int state_dim, int th );
void writeTimingFile(int run_idx, int world_size, long N, double elapsed,  double combining_accum, double rsamp_accum, double rsamp_accum2, short intact, int state_dim, int th);
void writeHistogram(double *x, long N, int run_idx);

#endif /* FILEIO_H_ */
