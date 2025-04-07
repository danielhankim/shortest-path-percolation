CC=gcc
OT=-ffast-math
OP=$(OT)
LOP=-lm -o

CHDR=./
CSRC=./
MSRC=./

MAIN = efficient_SPP_ERN
TAG=$(MAIN).out
OBJ=$(MAIN).o mt19937-64.o basic_functions.o percolation_functions.o file_io.o

$(TAG)	: $(OBJ)
	$(CC) $(OBJ) $(LOP) $(TAG)
	rm *.o

$(MAIN).o :
	$(CC) $(OP) -I$(CHDR) -c $(MSRC)$(MAIN).c

mt19937-64.o	:
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)mt19937-64.c

basic_functions.o	:
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)basic_functions.c

percolation_functions.o	:
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)percolation_functions.c

file_io.o	:
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)file_io.c

clear	:
	rm data/*; rm index.dat; cp work/words.lst ./

clean	:
	rm *.o
