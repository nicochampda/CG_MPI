CC=mpicc
CFLAGS=-O3
LIBS= -lmpi

conjgrad: conjgrad.c
	$(CC) $(CFLAGS) -o conjgrad conjgrad.c $(LIBS)

clean:
	$(RM) conjgrad
	$(RM) core.*
	$(RM) vgcore.*
