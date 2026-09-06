/* -------------------------------------------------------------------------- */
#include "eigen_serial.h"
/* -------------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


/* -------------------------------------------------------------------------- */
  

int main(int argc,char** argv) {
    // Generating random matrix, and printing it
    // srand(time(NULL));
    srand(10000);
    printf("Code started: \n");
    
    double qr_tol=1e-5;
    double max_iter=1000000;
    int n;
    if (argc>=2)
        n=atoi(argv[1]);
    else
    {
        n=400;
    }
    double *A = malloc(n*n*sizeof(double));
    double *Z = malloc(n*n*sizeof(double));
    double *d = malloc(n*sizeof(double));
    double *e = malloc(n*sizeof(double));
    char time_filename[100];
    sprintf(time_filename,"../data/time_n=%d_p=1.csv",n);
    FILE* f_time=fopen(time_filename,"w");
    fprintf(f_time,"n,p,k, Householder time, QL time, Total time\n");
    for (int k=0;k<10;k++)
    {   
        printf("k = %d\n",k);
        printf("n*n %d\n",n*n);
        char filename[100];


        sprintf(filename,"../data/initial_matrix_n=%d_%d.h5",n,k);
        printf("Filename = %s\n",filename);
        
        build_random_sym_matrix(A,n,-1,1);
        if (n<20)
        {    
            printf("A = \n");
            print_matrix(A,n);
        }
        // Saving initial matrix
        save_vec_hdf5(filename,A,n*n);
        
        time_t House_start = clock();
        householder(A,Z,d,e,n);
        time_t House_end = clock();

        time_t QL_start = clock();
        QL_iteration(Z,d,e,n,qr_tol,max_iter);
        time_t QL_end = clock();

        if (n<20)
        {    
            printf("Z = \n");
            print_matrix(Z,n);
            printf("d = \n");
            print_vector(d,n);
        }
        
        double duration=1.*(QL_end-House_start)/CLOCKS_PER_SEC;
        double House_duration=1.*(House_end-House_start)/CLOCKS_PER_SEC;
        double QL_duration=1.*(QL_end-QL_start)/CLOCKS_PER_SEC;
        printf("Time for one diagonalization = %4.4f s\n",duration);
        printf("Householder time = %4.4f s\n",House_duration);
        printf("QL time = %4.4f s\n",QL_duration);
        
        sprintf(filename,"../data/eigenvalues_n=%d_%d.h5",n,k);
        save_vec_hdf5(filename,d,n);
        sprintf(filename,"../data/eigenvectors_n=%d_%d.h5",n,k);
        save_vec_hdf5(filename,Z,n*n);
        // n,p,k, Householder time, QL time, Total time
        fprintf(f_time,"%d, %d, %d, %2.6f, %2.6f, %2.6f\n",n,1,k,House_duration,QL_duration,duration);
    }
    free(A);
    free(d);
    free(e);
    return 0;
}
