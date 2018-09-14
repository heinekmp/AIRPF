#include <mpi.h>

int 
findPairedProcessRank( 
		      int id, 
		      int s );

void createCommunicationGroups( 
			       int world_size, 
			       int S, 
			       MPI_Comm *group_comms, 
			       MPI_Group *groups, 
			       int world_rank );
void 
constructHypercube(
		  int world_size,
		  MPI_Comm * nthCube,
		  int * hc_rank,
		  int * nDim
		  );
void hypercube_matrix(
		      MPI_Comm * nthCube,
		      int *cube_mtx,
		      int S,
		      int world_rank
		      );

