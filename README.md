# AIRPF: Parallelising particle filters with butterfly interactions

The code included in this repository is not designed for general purposes and is merely provided for the transparency of the published material associated with this code. 

The success of running the code requires the directory structure to remain unchagned. Executables of the algorithms IPF1, IPF2, AIRPF1 and AIRPF2 can be created with the makefile with commands

```
make ipf1
make ipf2
make airpf1
make airpf2
```

after which the executables can be run with commands

```mpirun -np $NP ./ipf1 $I $N $D $TH
mpirun -np $NP ./ipf2 $I $N $D $TH
mpirun -np $NP ./airpf1 $I $N $D $TH
mpirun -np $NP ./airpf2 $I $N $D $TH
```
 
where the variables are as follows:

```
$NP   number of processors
$I    index for the run (any integer)
$N    number of particles to be used
$D    dimension of the state space
$TH   resampling threshold
```
Here a sample data is included in `data` folder which is of dimension 7. This means that one must use the parameter value `$D = 7` unless a new data with different dimension is provided.
