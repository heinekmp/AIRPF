#include <stdlib.h>

static long *shortage;
static long *excess;

void initialise_balancing(
			  int world_size) 
{
  shortage  = malloc( sizeof( long ) * world_size );
  excess = malloc( sizeof( long ) * world_size );
}

void finalise_balancing() 
{
  free(shortage);
  free(excess);
}

void 
calculate_balancing_matrix( 
			   long *cnts, 
			   int world_size, 
			   long *BM , 
			   long n ) 
{

  long d, s;
  
  for ( int i = 0; i < world_size ; i++ ) {
    s = n - cnts[ i ];
    shortage[ i ] = s > 0 ? s : 0;
  }
  
  for ( int i = 0; i < world_size ; i++ ) {
    d = cnts[ i ] - n;
    excess[ i ] = d > 0 ? d : 0;
  }

  for (int i = 0; i < world_size; i++ ) {
    for (int j = 0; j < world_size; j++ ) {
      
      if ( i == j ) {
	
	BM[ world_size * i + j ] = cnts[ i ];
	
      } else {
	
	if ( cnts[ j ] >= n ) { // excess 
	  
	  BM[ world_size * i + j ] = 0;
	  
	} else if ( shortage[ j ] > 0 ) { // shortage
	  
	  s = shortage[ j ]; // absolute shortage
	  
	  d = excess[ i ] < s ? excess[ i ] : s;
	  
	  excess[ i ] -= d;
	  
	  BM[ world_size * i + j ] = d;
	  
	  shortage[ j ] -= d;
	  
	} else {

	  BM[ world_size * i + j ] = 0;
	  
	}
      }   
    }
  }
}
