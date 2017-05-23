#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <math.h>
//#include <float.h>
#include <sys/time.h>
#include <unistd.h>

void initialize(double *A, double *b, int mat_size);
void matVectMult(double *A, double *p, int hei, int len, double *result);
double vectDotVect(double *v1, double *v2, int len);
void vectSum(double *u, double c, double *v, int len, double *result);
void Prvalues(int len, int hei,  double matrix[len * hei]);

int main(int argc, char *argv[]){
    int nprocs, rank;
    int nIter = 10;
    double max_error = 0.00001;
    int mat_size;

    double *A, *x;
    double *A_rows, *z, *r, *p, *p_loc;
    double alpha, gamma, gamma_new, gamma_loc, beta, p_z_loc;

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
    p = (double *)malloc(mat_size * sizeof(double));

    if(rank == 0){
        /* Initialization of A, b, and x */
        A = (double *)malloc(mat_size * mat_size * sizeof(double));
        initialize(A, p, mat_size);
        
        Prvalues(mat_size, mat_size, A);
        Prvalues(1, mat_size, p);


        /*Scatter input datas over MPI_COMM_WORLD */
        
        for (int k=0; k<nprocs; k++){
            MPI_Isend(&A[k*block_size*mat_size], block_size*mat_size, MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &init_send[k]);
        }
        MPI_Waitall(nprocs, init_send, MPI_STATUSES_IGNORE);
    }

    
    /* Last part of initialization */
    MPI_Bcast(p, mat_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    A_rows = (double *)malloc(block_size*mat_size*sizeof(double));

    MPI_Irecv(A_rows, block_size*mat_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &init_recv);

    r = (double *)malloc(block_size*sizeof(double));
    x = (double *)malloc(block_size*sizeof(double));
    z = (double *)malloc(block_size*sizeof(double));
    p_loc = (double *)malloc(block_size*sizeof(double));

    memset(x, 0, block_size * sizeof(double));
    memcpy(r, &p[rank * block_size], block_size*sizeof(double));

    MPI_Wait(&init_recv, MPI_STATUS_IGNORE);
    
    gamma = vectDotVect(p, p, mat_size);

    for(int n = 0; n<nIter; n++){
        /* Compute z=A*p */
        matVectMult(A_rows, p, block_size, mat_size, z);

        /* compute alpha and allreduce */
        p_z_loc = vectDotVect(&p[rank * block_size], z, block_size);
        MPI_Allreduce(&p_z_loc, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        alpha = gamma / alpha;
        
        /* compute x, r */
        vectSum(x, alpha, &p[rank * block_size], block_size, x);
        vectSum(r, -alpha, z, block_size, r);
        
        /* compute gamma and all reduce */
        gamma_loc = vectDotVect(r, r, block_size);
        MPI_Allreduce(&gamma_loc, &gamma_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        /* check for exiting loop */
        if (gamma_new < max_error){
            break;
        }
        //printf("gamma %f\n", gamma_new);
        /* compute beta in each procs */
        beta = gamma_new / gamma;
        gamma = gamma_new;

        /* compute p and allgather */
        vectSum(r, beta, &p[rank*block_size], block_size, p_loc);
        MPI_Allgather(p_loc, block_size, MPI_DOUBLE, p, block_size, MPI_DOUBLE, MPI_COMM_WORLD);


    } 
    
    sleep(rank);
    printf("rank %i\n", rank);
    //Prvalues(1, block_size, r);
    //Prvalues(mat_size, block_size, A_rows);
    //Prvalues(1, block_size, z);
    //printf("alpha : %f\n", alpha);
    Prvalues(1, block_size, x);
    //printf("gamma : %f\n", gamma_new);
    //Prvalues(1, block_size, p_loc);
    //Prvalues(1, mat_size, p);

    MPI_Finalize();

    return 0;
}

void initialize(double *A, double *b, int mat_size){
    double A_loc[mat_size*mat_size];
    for (int i=0; i<mat_size; i++){
        b[i] = rand()/(double) RAND_MAX;
        for (int j=0; j<mat_size; j++){
            A_loc[i*mat_size + j] = rand()/(double) RAND_MAX;
        }
        
    }
    for (int i=0; i<mat_size; i++){
        for(int j=0; j<mat_size; j++){
            for(int k=0; k<mat_size;k++){
                A[i*mat_size + j] += A_loc[i*mat_size + k] * A_loc[j*mat_size + k];
            }
        }
    }


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
    
void vectSum(double *u, double c, double *v, int len, double *result){
    for (int i=0; i<len; i++){
        result[i] = u[i] + c * v[i];
    }
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
