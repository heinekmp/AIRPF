/*
 * signal_data_generator.c
 *
 *  Created on: 5.5.2017
 *      Author: heine
 */

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "ziggurat.h"
#include "typedefs.h"

struct Process stochastic_volatility(long n) {

	float fn[128];
	uint32_t kn[128];
	float wn[128];
	uint32_t seed = 1;
	double phi = .9, w, sigma = 1;
	double *x, *y;
	struct Process process;

	r4_nor_setup(kn, fn, wn);

	x = malloc(sizeof(double) * n);
	y = malloc(sizeof(double) * n);

	x[0] = 0;
	y[0] =  exp(x[0]/2)*r4_nor(&seed, kn, fn, wn);

	for (long i = 1; i < n; i++) {
		w = r4_nor(&seed, kn, fn, wn);
		x[i] = phi * x[i - 1] + sigma*w;
		y[i] = exp(x[i]/2)*r4_nor(&seed, kn, fn, wn);
	}

	process.x = x;
	process.y = y;
	process.n = n;

	return process;

}
