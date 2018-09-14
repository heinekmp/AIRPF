#include <typedefs.h>

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
		 double *my_normal_weight,
		 int *do_resampling,
		 double rs_threshold
		 );
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
	     );

void 
compute_required_number_of_steps(
				 int world_size,
				 double rs_threshold,
				 int *N_steps,
				 double *proc_w,
				 int *cube_mtx,
				 int last_stage
				 );

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
		     );
