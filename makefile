IDIR = ./include
CC=mpicc
CFLAGS=-I$(IDIR) 
 
LIBS= -lm -lgsl -lgslcblas
ODIR = obj

_DEPS = fileio.h likelihood.h resampling.h typedefs.h mutation.h radix.h master_tasks.h randomisation.h ipf_balance.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ =  fileio.o likelihood.o resampling.o typedefs.o mutation.o radix.o master_tasks.o randomisation.o ipf_balance.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# All object files depend on the corresponding c-file and the DEPS
$(ODIR)/%.o: ./src/%.c $(DEPS)
	$(CC) -std=c99 -Wall -c -o $@ $< $(CFLAGS) $(LIBS)

# All the o-files need to be up-to-date
ipf2: ./obj/ipf2.o $(OBJ) 
	$(CC) -std=c99 -Wall -o $@ $^ $(CFLAGS) $(LIBS) 

# All the o-files need to be up-to-date
airpf1: ./obj/airpf1.o $(OBJ)
	$(CC) -std=c99 -Wall -o $@ $^ $(CFLAGS) $(LIBS)

ipf1: ./obj/ipf1.o $(OBJ)
	$(CC) -std=c99 -Wall -o $@ $^ $(CFLAGS) $(LIBS)

airpf2: ./obj/airpf2.o $(OBJ)
	$(CC) -std=c99 -Wall -o $@ $^ $(CFLAGS) $(LIBS)
