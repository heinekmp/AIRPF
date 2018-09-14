#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

static const int MASTER = 0;

int findPairedProcessRank(
			  int id, 
			  int s )   
{

  double pow20 = ( int ) pow( 2, s );
  double pow21 = ( double ) pow( 2, s + 1 );
  int pid1 = id % ( int ) pow20 + pow21 * floor( ( double ) id / pow21 );
  int pid2 = pid1 +  pow20;

  return pid1 == id ? pid2 : pid1;
 
}

/**********************************************/
/* HYPERCUBE CONSTRUCTION                     */
/**********************************************/
void 
constructHypercube(
		  int world_size,
		  MPI_Comm * nthCube,
		  int * hc_rank,
		  int * nDim
		  )
{		  
  
  // Dimension is the base to log of the world size
  * nDim = ( int ) log2( world_size );
  
  // Constant 2 elements per dimension
  int *processPerDim;
  processPerDim = malloc( sizeof( int ) * * nDim );
  for ( int i = 0; i < * nDim ; i++ )
    processPerDim[ i ] = 2;
  
  // Each dimension is periodical
  int *period;
  period = malloc( sizeof( int ) * * nDim );
  for ( int i = 0; i < * nDim ; i++ ) 
    period[ i ] = 1;
  
  // Create the hyper cube
  MPI_Cart_create( MPI_COMM_WORLD, * nDim, processPerDim, period, 1, nthCube );
  
  // My rank in the hypercube
  MPI_Comm_rank( * nthCube, hc_rank );
  
}

void hypercube_matrix(
		      MPI_Comm * nthCube,
		      int *cube_mtx,
		      int S,
		      int world_rank
		      )
{

  int src; // dummy
  int *pairs;

  pairs = malloc( sizeof( int ) * S );

  for ( int s = 0; s < S; s++) 

    MPI_Cart_shift( *nthCube, s, 1, &src, &pairs[s] );

  MPI_Gather( pairs, S, MPI_INT, cube_mtx, S, MPI_INT, MASTER, MPI_COMM_WORLD );

  free( pairs );
}
