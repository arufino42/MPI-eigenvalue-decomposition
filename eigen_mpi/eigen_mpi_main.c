/* -------------------------------------------------------------------------- */
#include "eigen_mpi.h"
/* -------------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <hdf5.h>
/* -------------------------------------------------------------------------- */

int main(int argc,char** argv) {
    // Generating random matrix, and printing it
    srand(10000);
    
    MPI_Init(&argc, &argv);
    int N_proc,rank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &N_proc);
    MPI_Comm_rank(comm, &rank);
    printf("Rank %d of %d\n",rank,N_proc);
    int n;
    if (argc>=2)
        n=atoi(argv[1]);
    else
    {
        n=400;
    }

    if (n%N_proc !=0 )
    {    
        printf("ERROR: Size of matrix should be divisible by the number of processors.\n");
        return 1;
    }

    
    int m=n/N_proc;    

    double qr_tol=1e-5;
    double max_iter=1000000;
    
    double *A = malloc(n*n*sizeof(double));
    double *Z = malloc(n*n*sizeof(double));
    double *Z_local = malloc(n*m*sizeof(double));
    double *A_local = malloc(n*m*sizeof(double));
    double *d = malloc(n*sizeof(double));
    double *e = malloc(n*sizeof(double));
    char time_filename[100];
    sprintf(time_filename,"../data_weak_scaling/time_n=%d_p=%d.csv",n,N_proc);
    FILE* f_time=NULL;
    if (rank==0)
    {
        f_time=fopen(time_filename,"w");
        fprintf(f_time,"n,p,k, Householder time, QL time, Total time\n");
    }
    for (int k=0;k<10;k++)
    {   
        printf("k = %d\n",k);
        printf("n = %d\n",n);
        char filename[100];
        sprintf(filename,"../data_weak_scaling/initial_matrix_mpi_n=%d_%d.h5",n,k);
        printf("Filename = %s\n",filename);
        
        if (rank==0)
        {
            // load_vec_hdf5(filename,A); 
            build_random_sym_matrix(A,n,-1,1);
            if (n<20)
            {    
                print_matrix(A,n,n);
            }
            time_t save_start = clock();
            save_vec_hdf5(filename,A,n*n);
            time_t save_end = clock();
            double save_duration=1.*(save_end-save_start)/CLOCKS_PER_SEC;
            printf("HDF5 save time = %4.4f s\n",save_duration);
        }
        // int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
        MPI_Scatter(A, n*m, MPI_DOUBLE,
            A_local, n*m, MPI_DOUBLE, 0,MPI_COMM_WORLD); 

        // Calling Householder method, performed in parallel in each local slice of the A matrix 
        time_t House_start = clock();
        householder(A_local,Z_local,d,e,m,n,rank);
        time_t House_end = clock();
        
        if (rank==0 && n<20)
        {
            printf("Z (Householder) = \n");
            print_matrix(Z_local,m,n);
        }
        
        time_t QL_start = clock();
        QL_iteration(Z_local,d,e,m,n,qr_tol,max_iter);
        time_t QL_end = clock();
        MPI_Allgather(Z_local,m*n,MPI_DOUBLE,Z,m*n,MPI_DOUBLE,MPI_COMM_WORLD);
        if (rank==0 && n<20)
        {
            printf("Z (end) = \n");
            print_matrix(Z,n,n);
            printf("Z_local (end) = \n");
            print_matrix(Z_local,m,n);
            printf("d = \n");
            print_vector(d,n);
        }
        if (rank==0)
        {
            double duration=1.*(QL_end-House_start)/CLOCKS_PER_SEC;
            double House_duration=1.*(House_end-House_start)/CLOCKS_PER_SEC;
            double QL_duration=1.*(QL_end-QL_start)/CLOCKS_PER_SEC;
            printf("Time for one diagonalization = %4.4f s\n",duration);
            printf("Householder time = %4.4f s\n",House_duration);
            printf("QL time = %4.4f s\n",QL_duration);
            sprintf(filename,"../data_weak_scaling/eigenvalues_mpi_n=%d_%d.h5",n,k);
            save_vec_hdf5(filename,d,n);
            sprintf(filename,"../data_weak_scaling/eigenvectors_mpi_n=%d_%d.h5",n,k);
            save_vec_hdf5(filename,Z,n*n);
            fprintf(f_time,"%d, %d, %d, %2.6f, %2.6f, %2.6f\n",n,N_proc,k,House_duration,QL_duration,duration);
        }
    }
    // Freeing memory and finalizing MPI
    free(A);
    free(A_local);
    free(Z);
    free(Z_local);
    free(d);
    free(e);
    MPI_Finalize();
    return 0;
}
