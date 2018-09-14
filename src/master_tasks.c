#define _XOPEN_SOURCE 600

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "fileio.h"

static const int MASTER = 0;

void 
combineEstimates(
		 int world_rank,
		 int state_dim,
		 struct Estimates estimates,
		 double * output,
		 long iteration,
		 int output_dim,
		 int world_size,
		 double * total_w,
		 double * my_normal_weight,
		 int * do_resampling,
		 double rs_threshold
)
{
     
  double w = 0.0;
  double w2;
  double new_w;
  double mean_w;
  double epc; // effective process count
  int do_rs = 1;
  
  double * new_x;
  new_x = malloc( sizeof( double ) * state_dim );

  double * est_buff;
  est_buff = malloc( sizeof( double ) * ( state_dim + 1 ) );

  int ests_offset = output_dim - world_size * state_dim;

  double * proc_w; // filter specific weights 
  proc_w = malloc( sizeof( double ) * world_size );
  
  double * estimate;
  estimate = malloc( sizeof( double ) * state_dim );

  if( world_rank > MASTER ) {
    // Process 0 is the combiner, others need to send their results
      
    // Send estimate to the process 0
    for ( int d = 0; d < state_dim; d++ )
      est_buff[ d ] = estimates.mu[ d ];

    est_buff[ state_dim ] = estimates.w;
    
    MPI_Send( est_buff, state_dim + 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
    
  } else {
    // Process 0 receives and combines the estimates
      
    for ( int d = 0; d < state_dim; d++ ) {
	
      estimate[ d ] = estimates.mu[ d ];
	
      output[ output_dim * iteration + ests_offset + d ] = estimate[ d ];
 
    }
    
    w = estimates.w;
    
    w2 = w * w;
      
    proc_w[ 0 ] = w;
    
    for( int i = 1; i < world_size; i++ ) {
      
      // Receive estimates with weights from other processes
      MPI_Recv( est_buff, state_dim + 1, MPI_DOUBLE, i, 0, 
		MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      
      for ( int d = 0; d < state_dim; d++ ) {
	
	new_x[ d ] = est_buff[ d ];
	
	output[ output_dim * iteration + ests_offset + i * state_dim + d ] = new_x[ d ];
	  
      }
	
      new_w = est_buff[ state_dim ];
	
      proc_w[ i ] = new_w;
	
      // combine estimates
      for ( int d = 0; d < state_dim; d++ )
	
	estimate[ d ] = ( estimate[ d ] * w + new_w * new_x[ d ] )  / ( w + new_w );
      
      w += new_w;
      
      w2 += new_w * new_w;
            
    }

    // Normalise process weights
    for (int i = 0; i < world_size; i++ )
      
      proc_w[i] = proc_w[i] / w;

    // Compute the output
    
    for ( int d = 0; d < state_dim; d++ )      
      output[ output_dim * iteration + d ] = estimate[ d ];
    
    output[ output_dim * iteration + state_dim ] = w;
    
    mean_w = w / ( double ) world_size;

    // normalised effective process count. Takes values in [1/world_size,1].
    epc = mean_w * mean_w  / ( w2 / ( double ) world_size );
    output[ output_dim * iteration + state_dim + 1 ] = epc;
    
    for (int i = 0; i < world_size; i++ ) 
      output[ output_dim * iteration + state_dim + 2 + i ] = proc_w[ i ];
    
    
    do_rs = epc < rs_threshold ? 1 : 0 ;  

    //    output[ output_dim * iteration + output_dim - 1 ] = do_rs * nDim;

  }
  
  // Send the total weight to all processes
  * total_w = w;
  MPI_Bcast( total_w, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
  
  
  // Broadcast resampling indicator to all processes
  MPI_Bcast( &do_rs, 1, MPI_INT, MASTER, MPI_COMM_WORLD );
  
  *do_resampling = do_rs;

  // Send each process their normalised weight 
  MPI_Scatter( proc_w, 1, MPI_DOUBLE, my_normal_weight, 1, MPI_DOUBLE, 
	       MASTER, MPI_COMM_WORLD );

  free( estimate );
  free( new_x );
  free( proc_w );
  free( est_buff );

}

// ----------------------
// Write the output files
// ----------------------
void 
createOutput(
	     int algorithm,
	     int world_rank,
	     double * output,
	     long data_length,
	     int world_size,
	     int run_idx,
	     long N,
	     int state_dim,
	     struct timespec start_time,
	     struct timespec end_time,
	     double comm_time,
	     double rsamp_time,
	     double comm_time2,
	     double weight_time,
	     double mutation_time,
	     double init_time,
	     int th
)
{
 
  if ( world_rank == MASTER ){
    
    writeOutputFile( output, data_length, world_size, run_idx, N, algorithm, state_dim,
		     th );
    
    double elapsed = ( ( double ) ( end_time.tv_sec - start_time.tv_sec ) * 1.0e9 + 
		      ( double ) ( end_time.tv_nsec - start_time.tv_nsec ) ) / 1.0e9;
    
    printf( "%10i %10i %10ld %10f %10f %10f %10f %10f %10f %10f %10i %10i\n", 
	    run_idx, world_size, N, elapsed, comm_time, rsamp_time, 
	    comm_time2, weight_time, mutation_time, init_time, algorithm, state_dim );
    
    writeTimingFile(run_idx, world_size, N, elapsed, 0, 0, 0, algorithm, state_dim, th);

  }
 
}

void 
compute_required_number_of_steps(
				 int world_size,
				 double rs_threshold,
				 int *N_steps,
				 double *proc_w,
				 int *cube_mtx,
				 int last_stage
				 )
{
  
  double mean_w = 0;
  double w2 = 0;
  int M = world_size;
  int S = ( int ) log2( world_size );
  double tmp;
  
  *N_steps = 0;
  
  // iterate over all possible step numbers (until sufficient number of steps 
  // is found)
  for ( int s = 0; s < S ; s++ ) {
    
    mean_w = 0;
    w2 = 0;

    // Calculate the effective process count
    for( int i = 0; i < M; i++ ) {
      
      mean_w += proc_w[ i ];
      w2 += proc_w[ i ] * proc_w[ i ];

    }
    
    mean_w = mean_w / M;

    // Compute and compare effective process count
    if (  mean_w * mean_w  / ( w2 / M )  < rs_threshold ) {
      
      *N_steps = *N_steps + 1;
      
      // Recompute the process weights
      for ( int j = 0; j < M; j++ ) {
	
	tmp = ( proc_w[ j ] + proc_w[ cube_mtx[ j * S + ( s + last_stage + 1 ) % S ] ] ) 
	  / (double) 2.0;
	
	proc_w[ j ] = tmp;
	proc_w[ cube_mtx[ j * S + ( s + last_stage + 1 ) % S ] ] = tmp;

      }
      
    } else {
      
      break;
      
    }
  }
  
}

void 
flexCombineEstimates(
		     int world_rank,
		     int state_dim,
		     struct Estimates estimates,
		     double * output,
		     long iteration,
		     int output_dim,
		     int world_size,
		     double * total_w,
		     double * my_normal_weight,
		     int *N_steps,
		     double rs_threshold,
		     int *cube_mtx,
		     int last_stage
		     )
{
     
  double w = 0.0;
  double w2;
  double new_w;
  double mean_w;
  double epc; // effective process count
  
  double * new_x;
  new_x = malloc( sizeof( double ) * state_dim );

  double * est_buff;
  est_buff = malloc( sizeof( double ) * ( state_dim + 1 ) );

  int ests_offset = output_dim - world_size * state_dim;

  double * proc_w; // filter specific weights 
  proc_w = malloc( sizeof( double ) * world_size );
  double *proc_w_flex;
  proc_w_flex = malloc( sizeof( double ) * world_size );
  
  double * estimate;
  estimate = malloc( sizeof( double ) * state_dim );

  if( world_rank > MASTER ) {
    // Process 0 is the combiner, others need to send their results
      
    // Send estimate to the process 0
    for ( int d = 0; d < state_dim; d++ )
      est_buff[ d ] = estimates.mu[ d ];

    est_buff[ state_dim ] = estimates.w;
    
    MPI_Send( est_buff, state_dim + 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
    
  } else {
    // Process 0 receives and combines the estimates
      
    for ( int d = 0; d < state_dim; d++ ) {
	
      estimate[ d ] = estimates.mu[ d ];
	
      output[ output_dim * iteration + ests_offset + d ] = estimate[ d ];
 
    }
    
    w = estimates.w;
    
    w2 = w * w;
      
    proc_w[ 0 ] = w;
    
    for( int i = 1; i < world_size; i++ ) {
      
      // Receive estimates with weights from other processes
      MPI_Recv( est_buff, state_dim + 1, MPI_DOUBLE, i, 0, 
		MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      
      for ( int d = 0; d < state_dim; d++ ) {
	
	new_x[ d ] = est_buff[ d ];
	
	output[ output_dim * iteration + ests_offset + i * state_dim + d ] = new_x[ d ];
	  
      }
	
      new_w = est_buff[ state_dim ];
	
      proc_w[ i ] = new_w;
	
      // combine estimates
      for ( int d = 0; d < state_dim; d++ )
	
	estimate[ d ] = ( estimate[ d ] * w + new_w * new_x[ d ] )  / ( w + new_w );
      
      w += new_w;
      
      w2 += new_w * new_w;
            
    }

    // Normalise process weights
    for (int i = 0; i < world_size; i++ ) {
      
      proc_w[ i ] = proc_w[ i ] / w;
      proc_w_flex[ i ] = proc_w[ i ];

    }

    compute_required_number_of_steps( world_size, rs_threshold, 
				      N_steps, proc_w_flex, cube_mtx, last_stage );
    
    output[ output_dim * iteration + output_dim - 1 ] = *N_steps;

    // Compute the output
    
    for ( int d = 0; d < state_dim; d++ )      
      output[ output_dim * iteration + d ] = estimate[ d ];
    
    output[ output_dim * iteration + state_dim ] = w;
    
    mean_w = w / ( double ) world_size;

    // normalised effective process count. Takes values in [1/world_size,1].
    epc = mean_w * mean_w  / ( w2 / ( double ) world_size );
    output[ output_dim * iteration + state_dim + 1 ] = epc;
    
    for (int i = 0; i < world_size; i++ ) 
      output[ output_dim * iteration + state_dim + 2 + i ] = proc_w[ i ];
    
  }
  
  // Send the total weight to all processes
  * total_w = w;
  MPI_Bcast( total_w, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD );
  
  // Broadcast required number of steps to all processes
  MPI_Bcast( N_steps, 1, MPI_INT, MASTER, MPI_COMM_WORLD );

  // Send each process their normalised weight 
  MPI_Scatter( proc_w, 1, MPI_DOUBLE, my_normal_weight, 1, MPI_DOUBLE, 
	       MASTER, MPI_COMM_WORLD );

  free( estimate );
  free( new_x );
  free( proc_w );
  free( est_buff );
  free( proc_w_flex );

}
