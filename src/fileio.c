/*
 * fileio.c
 *
 *  Created on: 17.5.2017
 *      Author: heine
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <float.h>
#include "resampling.h"
#include "typedefs.h"

const char *comb_est_file = "comb_est.csv";
const char *process_file_x = "process_x.csv";
const char *process_file_y = "./data/process_y_";
const char *data_extension = ".csv";
const char *result_folder = "results";
const char *folder = "fileexchange";
const char *outfile = "output.csv";
const char *processfile = "processes.csv";
const char *result_outfile = "estimate.csv";
long histogram[1000];

struct Process readProcessFile( int state_dim ) {
  
  FILE *file;
  char *token, line[ 1000 ];
  char filename[ 80 ];
  const char s[ 2 ] = ";";
  const long max_length = 10000;
  const int max_dim = 100;
  long it_count = 0;
  struct Process process;
  double *tmp;
  int dim;
  
  tmp = malloc( sizeof( double ) * max_length * max_dim );
  
  process.n = 0;
  
  it_count = 0; 
  
  sprintf( filename, "%s%i%s", process_file_y, state_dim, data_extension );
  
  file = fopen( filename, "r" );
  
  if ( file != NULL ) {
    
    while( fgets( line, sizeof line, file ) != NULL ) {
      
      token = strtok(line, s);
      
      dim = 0;
      
      // iterate over dimensions
      while ( token != NULL ) {
	
	tmp[ max_length * dim + it_count ] = atof(token);
	
	token = strtok(NULL, s);
	
	dim++;
	
      }
      
      if (it_count == 0) 
	
	process.dim = dim;
      
      it_count++;
      
    }	
    
  } else {
    
    printf("Unable to read the process file.\n");
    
  }
  
  process.n = it_count;
  
  fclose(file);
  
  process.y = malloc( sizeof( double ) * it_count * process.dim );
  
  for (int it = 0; it < it_count; it++ ) {
    for (int d = 0; d < process.dim; d++ ) {
      process.y[ it_count * d + it ] = tmp[ max_length * d + it ];
    }
  }
  
  free(tmp);
  
  return process;
}


void reset_output() {
  remove(outfile);
  remove(processfile);
}

/*void writeProcessFile(struct Process process, int id) {

  FILE *file;
  char filename[80];
  sprintf(filename, "%i_%s", id, processfile);
  file = fopen(filename, "w");
  
  for (long i = 0; i < process.n; i++)
    fprintf(file, "%.16f;%.16f\n", process.x[i], process.y[i]);
  
  fclose(file);

  }*/

void output(double *w, double *x, long N, int id, long k, int s) {
  
  FILE *file;
  char tmpfilename[80], filename[80];
  
  sprintf(filename, "./%s/%i_%ld_%i_%s", folder, id, k, s, outfile);
  sprintf(tmpfilename, "./%s/%i_tmp.csv", folder, id);
  file = fopen(tmpfilename, "w");
  
  for (long i = 0; i < N - 1; i++) {
    fprintf(file, "%.16e;%.16e\n", w[i], x[i]);
  }
  fprintf(file, "%.16e;%.16e\n", w[N - 1], x[N - 1]);
  fclose(file);
  
  rename(tmpfilename, filename);
}

long countlines(char *filename) {
  
  FILE *file;
  file = fopen(filename, "r");
  char buf[80];
  long count = 0;
  while (fgets(buf, sizeof(buf), file) != NULL) {
    count++;
  }
  return count;
}

void combineSamples(double *like, double *x, long N, int id, long k, int stage,
		    double *comb_like, double *comb_x) {
  
  FILE *file;
  char filename[80], *token, line[50];
  const char s[1] = ";";
  
  long j;
  
  sprintf(filename, "./%s/%i_%ld_%i_%s", folder, id, k, stage, outfile);
  file = fopen(filename, "r");
  
  // Fill the first half of the arrays
  for (j = 0; j < N; j++) {
    comb_like[j] = like[j];
    comb_x[j] = x[j];
  }
  
  j = N;
  
  // Fill the second half of the vectors
  for (j = N; j < 2 * N; j++) {
    
    fgets(line, sizeof line, file);
    token = strtok(line, s);
    
    for (int i = 0; i < 2; i++) {
      if (i == 0) {
	comb_like[j] = atof(token);
	token = strtok(NULL, s);
      } else {
	comb_x[j] = atof(token);
      }
    }
	}
  fclose(file);
}

void removeOldFiles(long k, int nid, int S) {
  char filename[80];
  for (int id = 0; id < nid; id++)
    for (int s = 0; s < S; s++) {
      sprintf(filename, "./%s/%i_%ld_%i_%s", folder, id, k - 1, s, outfile);
      remove(filename);
    }
}

struct Estimates estimateOutput(double *w, double *x, int id, long N, int run_idx, int state_dim) {

  double *estimate;
  double weight = 0;
  struct Estimates estimates;
  
  estimate = malloc( sizeof( double ) * state_dim );


  for ( int d = 0; d < state_dim; d++ )

    estimate[ d ] = 0;


  for ( long i = 0; i < N; i++) {
    
    for ( int d = 0; d < state_dim; d++ ) 

      estimate[ d ] += w[ i ] * x[ state_dim * i + d ];

    weight += w[i];

  }


  for ( int d = 0; d < state_dim; d++ ) 
  
    estimate[ d ] = estimate[ d ] / weight;
  

  estimates.nzer = weight;
  estimates.mu = estimate;
  estimates.w = weight;

  return estimates;

}

void read_true_normalisers(double *mu) {

	FILE *file;
	char line[50];

	file = fopen("true.txt", "r");

	for (int j = 0; j < 20; j++) {
		fgets(line, sizeof line, file);
		mu[j] = atof(line);
	}
}
void read_true_normaliser_vars(double *sd) {

	FILE *file;
	char line[50];

	file = fopen("true_var.txt", "r");

	for (int j = 0; j < 20; j++) {
		fgets(line, sizeof line, file);
		sd[j] = atof(line);
	}
}

void writeIndices(long *inds, double *w, long N) {
	FILE *file;
	file = fopen("inds.txt", "w");
	for (long i = 0; i < N; i++)
		fprintf(file, "%ld %.16f\n", inds[i], w[i]);
	fclose(file);
}

void combinedEstimateOutput(double x, double w, int run_idx, int world_size, long N) {
  
  FILE *file;
  char filename[80];

  sprintf(filename, "./%s/%i_%i_%.0e%s", result_folder, world_size,
	  run_idx, (double)N, comb_est_file);

  file = fopen(filename,"a");
  
  fprintf(file,"%.10e;%.10e\n",x,w);

  fclose(file);
  
}

void writeOutputFile( 
		     double *data, 
		     long data_length, 
		     int world_size, 
		     int run_idx, 
		     long N, 
		     short intact, 
		     int state_dim, 
		     int th 
		      ) 
{
  
  FILE *file;
  char filename[ 80 ];
  int output_dim = state_dim + 2 + 2 * world_size + world_size * state_dim ;
  
  sprintf( filename, "./%s/%i_%i_%.7e_%i_%i_%i%s", result_folder, world_size, 
	   run_idx, ( double ) N, ( int ) intact, state_dim, th, comb_est_file );
  
  file = fopen( filename, "a" );

  for ( long i = 0; i < data_length; i++ ) {

    // Total weight across processes and average weight per process
    fprintf( file, "%.10e;%.10e;", data[ output_dim * i  + state_dim ], 
	     data[ output_dim * i  + state_dim + 1] );

    // Combined estimate per state dim
    for( int d = 0; d < state_dim; d++ )      
      fprintf( file, "%.10e;", data[ output_dim * i + d ] );

    fprintf( file, "%.0f;", data[ output_dim * i + output_dim - 1] );

    // // Process weights
    // for( int j = 0; j < world_size; j++ )      
    //  fprintf( file, "%.10e;", data[ output_dim * i + state_dim + 2 + j] );

    // // Process-wise resampling counts
    // for( int j = 0; j < world_size; j++ )
    // fprintf( file, "%.0f;", data[ output_dim * i + state_dim + 2 + world_size + j ] );

    // // Process-wise state estimates per dimension
    // for( int j = 0; j < world_size; j++ )
    //  for ( int d = 0; d < state_dim; d++ )
    // fprintf( file, "%.2f;", 
    //		 data[ output_dim * i + state_dim + 2 + 2 * world_size + 
    //		       j * state_dim + d ] );
    
    fprintf( file, "\n" );
  }

  fclose( file );
}

void writeTimingFile(int run_idx, int world_size, long N, double elapsed,  double combining_accum, double rsamp_accum, double rsamp_accum2, short intact, int state_dim, int th) {
  
  FILE *file;
  char filename[80];
  
  sprintf(filename,"./%s/%i_%i_%.7e_%i_%i_%i%s", result_folder, world_size, run_idx, (double)N, (int)intact, state_dim, th, "timing.csv");

  file = fopen(filename,"w");

  fprintf(file, "%i\t%i\t%ld\t%f\t%f\t%f\t%f\n", run_idx, world_size, N, elapsed, 
	  combining_accum, rsamp_accum, rsamp_accum2);

  fclose(file);

}

void writeHistogram(double *x, long N, int run_idx) {
  
  FILE *file;
  char filename[80];
  int bin;

  // Reset the histogram
  for(int i = 0; i < 1000; i++ )
    histogram[i] = 0;

  // Compute the histogram
  for(long i = 0; i < N; i++){
    bin = (int)( ( x[i] + 6.0 ) * 83.25 );
    bin = bin < 0 ? 0 : bin > 999 ? 999 : bin;
    histogram[bin]++;
  }

  // Write the histogram in a file
  sprintf(filename,"./%s/%i_%.7e_histogram.csv", result_folder, run_idx, (double)N );

  file = fopen(filename,"a");
  
  for(int i = 0; i < 1000; i++) {
    fprintf(file,"%ld;",histogram[i]);
  }
  fprintf(file,"\n");
  fclose(file);

}
