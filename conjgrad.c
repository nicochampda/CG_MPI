#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <float.h>
#include <sys/time.h>
//#include <unistd.h>


int main(int argc, char *argv[]){
    int nprocs, rank;
    int nIter = 1;
    int length, heigth;

    double *A, *b, *x;
    double *z, *r, *p;
    double alpha, gamma, gamma_new, beta;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    

    if(rank == 0){
        /* Initialization of A, b, and x */
        A = (double *)malloc(heigth * length * sizeof(double));
        b = (double *)malloc(heigth * sizeof(double));
        
        /*Scatter input datas over MPI_COMM_WORLD */

    }

    
    /* Last part of initialization */


    for(int n = 0; n<nIter; n++){
        /* Compute z=A*p */

        /* compute alpha and allreduce */

        /* compute x, r */

        /* compute gamma and all reduce */

        /* check for exiting loop */

        /* compute beta in each procs */

        /* compute p and allgather */
    } 


    MPI_Finalize();

    return 0;
}
