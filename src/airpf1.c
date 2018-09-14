#define _XOPEN_SOURCE 600

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <errno.h>
#include "fileio.h"
#include "likelihood.h"
#include "master_tasks.h"
#include "mutation.h"
#include "radix.h"
#include "randomisation.h"
#include "resampling.h"
#include "typedefs.h"

static const int ALG = 6;
static const int MASTER = 0;

int 
main (
      int argc, 
      char *argv[]) 
{

  struct Process process;
  struct Estimates estimates;
  int pair_rank, src;
  int do_resampling;
  long *inds;
  double *w, *x, *x_resamp;
  double total_w = 1;
  double my_normal_weight;
  double w_pair[2];
  double rsamp_time=0;
  double rsamp_start=0;
  double comm_time=0;
  double comm_start;
  double comm_time2=0;
  double *output;
  double weight_time = 0;
  double mutation_time = 0.0;
  double init_time = 0.0;
   unsigned int cnts[2];
  short counts[2];
  short n_resamp;
   
  const gsl_rng_type * T;
  double timer_start;

  // Parse commandline arguments
  int run_idx = 0;
  if ( argc > 1 )
    run_idx = atoi( argv[ 1 ] );

  long N = 1000;
  if ( argc > 2 )
    N = atol( argv[ 2 ] );

  int state_dim;
  if ( argc > 3 ) 
    state_dim = atoi( argv[ 3 ] );
  else {
    perror("State dimesion must be provided");
    return -1;
  }

  double epc_threshold;
  if (argc > 4 ) {
    epc_threshold = atof( argv[ 4 ] );
  } else {
    perror("Reampling threshold must be provided");
    return -1;
  }

  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
 
  gsl_rng_env_setup();
  T = gsl_rng_default;
  gsl_rng *rgen;
  rgen = gsl_rng_alloc( T );
 
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
 
  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Start the timer
  struct timespec start_time;
  clock_gettime(CLOCK_REALTIME, &start_time);
        
  MPI_Barrier( MPI_COMM_WORLD );
  timer_start = MPI_Wtime();

  long data_length;
  process = readProcessFile( state_dim );
  data_length = (int)process.n;
    
  // Construct the communication hypercube
  MPI_Comm nthCube;
  int hc_rank;
  int nDim;
  constructHypercube( world_size, &nthCube, &hc_rank, &nDim );
  int S = nDim;

  //---------------
  // Initialisation
  //---------------
  x = malloc( sizeof( double ) * N * state_dim ); 
  w = malloc( sizeof( double ) * N );

  inds = malloc( sizeof( long ) * N );
  x_resamp = malloc( sizeof( double ) * N * state_dim );
   
  // Initial distribution setup
  setupInitialSeed( world_rank, rgen );

  // Create initial sample
  createInitialSample( x, N, rgen, state_dim );
  
  // Setup the randomisation
  set_resamplingseed( world_rank );

  int output_dim = 2 + state_dim + 2 * world_size + world_size * state_dim;
  if ( world_rank == MASTER )
    output = malloc( sizeof( double ) * output_dim * data_length );
  
  my_normal_weight = 1.0 / world_size;
 
  MPI_Barrier( MPI_COMM_WORLD );
  init_time += MPI_Wtime() - timer_start;

  // main loop
  for( long n = 0; n < data_length; n++ ) {

    // ------------------------------------------------------------------------------
    MPI_Barrier( MPI_COMM_WORLD );
    timer_start = MPI_Wtime();

    evaluateLikelihoods( x, process, state_dim, data_length, N, n, w, 
			 my_normal_weight ); 
    
    MPI_Barrier( MPI_COMM_WORLD );
    weight_time += MPI_Wtime() - timer_start ;

    // Compute the relevant estimates for THIS filter only
    estimates = estimateOutput( w, x, world_rank, N, run_idx, state_dim );
    
    // Combine the results of all filters
    combineEstimates( world_rank, state_dim, estimates, output, n, 
		      output_dim, world_size, &total_w, &my_normal_weight, 
		      &do_resampling, epc_threshold);
    
    // -----------------------------
    // RESAMPLING
    // -----------------------------

    // Resampling within processors

    MPI_Barrier( MPI_COMM_WORLD );
    rsamp_start = MPI_Wtime();
    
    serial_multinomial_resample_sng_noutk( N, w, inds, N );
    
    for ( long i = 0; i < N; i++)
      for ( int d = 0; d < state_dim; d++ ) 
	x_resamp[ state_dim * i + d ] = x[ state_dim * inds[ i ] + d ];

    for ( long i = 0; i < N; i++)
      for ( int d = 0; d < state_dim ; d++ ) 
	x[ state_dim * i + d ] = x_resamp[ state_dim * i + d ];
    
    // -----------------------------------------------------------------------------
    MPI_Barrier( MPI_COMM_WORLD );
    rsamp_time += MPI_Wtime() - rsamp_start;
    comm_start = MPI_Wtime();
    
    if ( do_resampling ) {

      // ---------------------------
      // START THE RADIX INTERACTION
      // ---------------------------
      
      for ( int s = 0; s < S ; s++ ) {
        
	// Determine the process you need to communicate with
	MPI_Cart_shift( nthCube, s, 1, &src, &pair_rank );
	
	// The lower ranking process in the pair gets the filter weights and
	// does the island level resampling 
	if ( hc_rank < pair_rank ) {
	  
	  MPI_Recv( w_pair + 1, 1, MPI_DOUBLE, pair_rank, 0, nthCube, 
		    MPI_STATUS_IGNORE );
	  
	  // The lower ranking filter must send its weight to the paired filter 
	  // so that the paired filter can update its weight
	  MPI_Send( &my_normal_weight, 1, MPI_DOUBLE, pair_rank, 1, nthCube ); 
	  
	  w_pair[ 0 ] = my_normal_weight;
	  
	  // Generate a binomial random variate indicating how many samples 
	  // are drawn from each of the paired processes.
	  gsl_ran_multinomial( rgen , ( size_t ) 2, ( unsigned int ) ( 2 ), w_pair, 
			       cnts );

	  counts[ 0 ] = ( short ) cnts[ 0 ];
	  counts[ 1 ] = ( short ) cnts[ 1 ]; 
	  
	  // Send the resampling counts to the paired process
	  MPI_Send( counts + 1 , 1, MPI_SHORT, pair_rank, 3, nthCube );
	  
	  n_resamp = counts[ 0 ];
	  
	} else { 
	  
	  // Send my weight to the lower ranking filter in the pair
	  MPI_Send(&my_normal_weight, 1, MPI_DOUBLE, pair_rank, 0, nthCube );
	  
	  // Receive the weight of the paired filter
	  MPI_Recv( w_pair, 1, MPI_DOUBLE, pair_rank, 1, nthCube, MPI_STATUS_IGNORE );
	  
	  // Receive the duplicate count ( 0 or 1 )
	  MPI_Recv( &n_resamp, 1, MPI_SHORT, pair_rank, 3, nthCube, MPI_STATUS_IGNORE );
	  
	  w_pair[ 1 ] = my_normal_weight;
	  
	}
	
	// -----------------------------------------------------------------------------  
	MPI_Barrier( MPI_COMM_WORLD );
	comm_time2 += MPI_Wtime() - comm_start;
	comm_start = MPI_Wtime();
	
	// Balance the resampled particle counts
	if ( n_resamp > 1 ) {
	  
	  // Send the surplus particles, i.e. the last n_resamp - N of the 
	  // resampled particles
	  MPI_Send( x, N * state_dim , MPI_DOUBLE, pair_rank, 2, nthCube );
	  
	} else if ( n_resamp < 1 ) { 
	  
	  // Receive compensation for shortage
	  MPI_Recv( x,  N * state_dim, MPI_DOUBLE, pair_rank, 2, nthCube, 
		    MPI_STATUS_IGNORE );
	  
	}
	
	my_normal_weight = ( w_pair[ 0 ] + w_pair[ 1 ] ) / (double ) 2;
      }

    }
      
    // -----------------------------------------------------------------------------
    MPI_Barrier( MPI_COMM_WORLD );
    comm_time += MPI_Wtime() - comm_start;
    timer_start = MPI_Wtime();

    mutation(N, x, state_dim, rgen);
    
    // -----------------------------------------------------------------------------
    MPI_Barrier( MPI_COMM_WORLD );
    mutation_time += MPI_Wtime() - timer_start ;

  }
  

  struct timespec end_time;
  clock_gettime(CLOCK_REALTIME, &end_time);
  
  createOutput( ALG, world_rank, output, data_length, world_size,
		run_idx, N, state_dim, start_time, end_time, comm_time, rsamp_time,
		comm_time2, weight_time,  mutation_time, init_time, 
		(int)( atof( argv[ 4 ] ) * 100 ) );
  
  // Finalize the MPI environment.
  MPI_Finalize();
  
  return 0;
}
