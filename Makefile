CC=mpicc
CFLAGS=-Wall -O3 -g
LIBS= -lmpi

conjgrad: conjgrad.c
	$(CC) $(CFLAGS) -o conjgrad conjgrad.c $(LIBS)

clean:
	$(RM) conjgrad
	$(RM) core.*
