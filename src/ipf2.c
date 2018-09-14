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
#include "ipf_balance.h"
#include "likelihood.h"
#include "master_tasks.h"
#include "mutation.h"
#include "randomisation.h"
#include "resampling.h"
#include "typedefs.h"

static const int ALG = 5;
static const int MASTER = 0;

int 
main (
      int argc, 
      char *argv[] ) 
{

  struct Process process;
  struct Estimates estimates;
  long * inds;
  long * BM;
  int do_resampling;
  double * w;
  double * x;
  double * x_resamp;
  double * proc_w;
  double * output;
  double my_normal_weight;
  double rsamp_time = 0;
  double rsamp_start = 0;
  double comm_time = 0;
  double comm_start;
  double comm_time2 = 0;
  double total_w = 1.0;
  double weight_time = 0.0;
  double mutation_time = 0.0;
  double init_time = 0.0;
  long n_xfer;
  unsigned int *cnts;
  long *cnts_l;
  const gsl_rng_type *T;
  double timer_start;

  // Parse command line arguments
  int run_idx = 0;
  if(argc > 1)
    run_idx = atoi(argv[1]);

  long N = 1000;
  if(argc > 2)
    N = atol(argv[2]);

  int state_dim;
  if ( argc > 3 ) 
    state_dim = atoi( argv[ 3 ] );
  else { 
    perror("State dimesion must be provided");
    return -1;
  }

  double epc_threshold;
  if (argc > 4 ) 
    epc_threshold = atof( argv[ 4 ] );
  else {
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
  struct timespec start_time;;
  clock_gettime(CLOCK_REALTIME, &start_time);

  MPI_Barrier( MPI_COMM_WORLD );
  timer_start = MPI_Wtime();
  
  long data_length;
  process = readProcessFile( state_dim );
  data_length = (int)process.n;
 
  //---------------
  // Initialisation
  //---------------
  x = malloc( sizeof( double ) * N * state_dim); 
  w = malloc( sizeof( double ) * N ); 

  inds = malloc( sizeof( long ) * N );
  x_resamp = malloc( sizeof( double ) * N * state_dim);
    
  proc_w = malloc( sizeof( double ) * world_size );
  cnts = malloc( sizeof( unsigned int ) * world_size );
  cnts_l = malloc( sizeof( long ) * world_size );

  BM = malloc( sizeof( long ) * world_size * world_size );
  initialise_balancing( world_size );

  // Initial distribution setup 
  setupInitialSeed( world_rank, rgen );

  // Create initial sample
  createInitialSample( x, N, rgen, state_dim );

  // Setup the randomisation
  set_resamplingseed( world_rank );

  int output_dim = 2 + state_dim + 2 * world_size + world_size * state_dim;
  if ( world_rank == MASTER )
    output = malloc( sizeof( double ) * output_dim * data_length );

  MPI_Barrier( MPI_COMM_WORLD );
  init_time += MPI_Wtime() - timer_start;

  my_normal_weight = 1.0 / world_size;

  // main loop
  for( long n = 0; n < data_length; n++ ) {

    // -----------------------------------------------------------------------------
    MPI_Barrier( MPI_COMM_WORLD );
    timer_start = MPI_Wtime();

    evaluateLikelihoods( x, process, state_dim, data_length, N, n, w, 
			 my_normal_weight );
    
    MPI_Barrier( MPI_COMM_WORLD );
    weight_time += MPI_Wtime() - timer_start ;

    // compute the relevant estimates
    estimates = estimateOutput( w, x, world_rank, N, run_idx, state_dim );
    
    // Combine the results of all filters
    combineEstimates( world_rank, state_dim, estimates, output, n,
		      output_dim, world_size, &total_w, &my_normal_weight,
		      &do_resampling, epc_threshold );
    //    printf("My weight %i: %f\n", (int)world_rank,my_normal_weight);
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
   
      //my_weight = estimates.w; // sum of weights for the current process

      // --------------------------
      // START THE FULL INTERACTION
      // --------------------------
      
      // Processes send their weights to the MASTER process
      MPI_Gather( &my_normal_weight, 1, MPI_DOUBLE, proc_w, 1, MPI_DOUBLE, MASTER, 
		  MPI_COMM_WORLD );
      MPI_Barrier( MPI_COMM_WORLD );
      
      // Sample how many particles need to be taken from each process
      if ( world_rank == MASTER ) {
	
	gsl_ran_multinomial( rgen, ( size_t ) world_size, 
			     ( unsigned int ) ( world_size ), proc_w, cnts );
	
	for( int j = 0; j < world_size; j++ )
	  cnts_l[ j ] = ( long ) cnts[ j ];
	
      }
      
      // Send the island resampled index array to all processes
      MPI_Bcast( cnts_l, world_size, MPI_LONG, MASTER, MPI_COMM_WORLD ); 
      
      //if ( world_rank == MASTER ) 
      //  for ( int j = 0; j < world_size; j++ ) 
      //    output[ output_dim * n + count_offset + j ] = cnts_l[ j ];
      
      // -----------------------------------------------------------------------------
      MPI_Barrier( MPI_COMM_WORLD );
      comm_time2 += MPI_Wtime() - comm_start;
      comm_start = MPI_Wtime();
      
      // Calculate how the islands should be communicated
      calculate_balancing_matrix( cnts_l, world_size , BM , 1 );
      
      if ( cnts_l[ world_rank ] < 1 ) { // this process needs an island
	
	// iterate over all processes (other than the current) and see if that 
	// process has surplus particles
	for ( int i = 0; i < world_size; i++ ) {
	  
	  if ( i != world_rank ) {
	    
	    // NB! BM[ world_size * i + world_rank ] should be 0 or 1
	    n_xfer = BM[ world_size * i + world_rank ] * N; 
	    
	    if ( n_xfer > 0 ) {
	      
	      MPI_Recv( x, N * state_dim, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, 
			MPI_STATUS_IGNORE );
	      
	      cnts_l[ world_rank ] += n_xfer;
	      
	    }
	  }
	}
	
      } else if ( cnts_l[ world_rank ] > 1 ) { // this process has excess islands
	
	// iterate over all processes (other than the current) and see if 
	// they need particles
	for (int i = 0; i < world_size; i++ ) {
	  if( world_rank != i ) {
	    
	    // NB! BM[ world_size * i + world_rank ] should be 0 or 1
	    n_xfer = BM[ world_size * world_rank + i] * N;
	    
	    if ( n_xfer > 0 ) {
	      
	      MPI_Send( x, n_xfer * state_dim, MPI_DOUBLE, i, 2, MPI_COMM_WORLD );
	      
	      cnts_l[ world_rank ] -= n_xfer;
	      
	    }
	  }
	}
	
      }

      my_normal_weight = 1.0 / world_size;
    } 
    //else {
    //      printf("Skip resample!\n");
    //    }

    //if ( world_rank == MASTER ) {
    //  printf("%i: %f\n", (int)n,my_normal_weight);
    //}


    // ----------------------------------------------------------------------------
    MPI_Barrier( MPI_COMM_WORLD );
    comm_time += MPI_Wtime() - comm_start;
    timer_start = MPI_Wtime();

    mutation(N, x, state_dim, rgen );

    // -----------------------------------------------------------------------------
    MPI_Barrier( MPI_COMM_WORLD );
    mutation_time += MPI_Wtime() - timer_start ;
    
  }  
  
  struct timespec end_time;
  clock_gettime(CLOCK_REALTIME, &end_time);
  
  createOutput( ALG, world_rank, output, data_length, world_size,
		run_idx, N, state_dim, start_time, end_time, comm_time, rsamp_time, 
		comm_time2, weight_time, mutation_time, init_time, 
		(int)( atof( argv[ 4 ] ) * 100 ) );

  finalise_balancing();

  // Finalize the MPI environment.
  MPI_Finalize();

  return 0;
}


