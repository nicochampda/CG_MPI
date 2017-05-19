#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <math.h>
//#include <float.h>
#include <sys/time.h>
#include <unistd.h>


void matVectMult(double *A, double *p, int hei, int len, double *result);
double vectDotVect(double *v1, double *v2, int len);
double *vectSum(double *u, double c, double *v, int len);
void Prvalues(int len, int hei,  double matrix[len * hei]);

int main(int argc, char *argv[]){
    int nprocs, rank;
    int nIter = 1;
    int mat_size;
    int i, j, k;

    double *A, *x;
    double *A_rows, *z, *r, *p;
    double alpha, gamma, gamma_new, beta, p_z_loc;

    if (argc != 2){
        printf("usage : mpirun -np nprocs ./conjgrad mat_size");
        return -1;
    }

    mat_size = atoi(argv[1]);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Request init_send[nprocs];
    MPI_Request init_recv;
    
    int block_size = mat_size/nprocs;
    r = (double *)malloc(mat_size * sizeof(double));

    if(rank == 0){
        /* Initialization of A, b, and x */
        A = (double *)malloc(mat_size * mat_size * sizeof(double));

        for (i=0; i<mat_size; i++){
            r[i] = i;
            for (j=0; j<mat_size; j++){
                A[i*mat_size + j] = i + j;
            }
        
        }
        /*Scatter input datas over MPI_COMM_WORLD */
        
        for (k=0; k<nprocs; k++){
            MPI_Isend(&A[k*block_size*mat_size], block_size*mat_size, MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &init_send[k]);
        }
        MPI_Waitall(nprocs, init_send, MPI_STATUSES_IGNORE);
    }

    
    /* Last part of initialization */
    MPI_Bcast(r, mat_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    A_rows = (double *)malloc(block_size*mat_size*sizeof(double));

    MPI_Irecv(A_rows, block_size*mat_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &init_recv);

    p = (double *)malloc(mat_size*sizeof(double));
    x = (double *)malloc(block_size*sizeof(double));
    z = (double *)malloc(block_size*sizeof(double));

    memset(x, 0, block_size * sizeof(double));
    memcpy(p, r, mat_size*sizeof(double));

    MPI_Wait(&init_recv, MPI_STATUS_IGNORE);
    
    gamma = vectDotVect(r, r, mat_size);

    for(int n = 0; n<nIter; n++){
        /* Compute z=A*p */
        matVectMult(A_rows, p, block_size, mat_size, z);

        /* compute alpha and allreduce */
        p_z_loc = vectDotVect(&p[rank * block_size], z, block_size);
        MPI_Allreduce(&p_z_loc, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        alpha = gamma / alpha;
        
        /* compute x, r */
        x = vectSum(x, alpha, &p[rank * block_size], block_size);
        r = vectSum(r, -alpha, z, block_size);
        /* compute gamma and all reduce */

        /* check for exiting loop */

        /* compute beta in each procs */

        /* compute p and allgather */
    } 
    
    sleep(rank);
    printf("rank %i\n", rank);
    Prvalues(1, mat_size, r);
    Prvalues(mat_size, block_size, A_rows);
    Prvalues(1, block_size, z);
    printf("alpha : %f\n", alpha);
    Prvalues(1, block_size, x);

    MPI_Finalize();

    return 0;
}


void matVectMult(double *mat, double *v, int hei, int len, double *result){
    for(int i=0; i<hei; i++){
        result[i] = vectDotVect(&mat[i*len], v, len);
    }
}

double vectDotVect(double *v1, double *v2, int len){
    double result = 0;
    for(int i=0; i<len; i++){
        result += (v1[i] * v2[i]);
    }
    return result;
}    
    
double *vectSum(double *u, double c, double *v, int len){
    for (int i=0; i<len; i++){
        u[i] = u[i] + c * v[i];
    }
    return u;
}


void Prvalues(int len, int hei,  double matrix[len * hei]){   
    int i, j;
    printf("\n");
    for (i = 0; i < hei; i++){
        for (j = 0; j < len; j++){
            printf("%.3f\t", matrix[i*len + j]);
        }
        printf("\n");
    }
    printf("\n");
}
