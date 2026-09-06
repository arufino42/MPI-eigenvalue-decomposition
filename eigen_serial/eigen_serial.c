#include "eigen_serial.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void save_vec_hdf5(const char *file, double *v, int len)
{
    hid_t f = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t d[1] = { len };
    hid_t s = H5Screate_simple(1, d, NULL);

    hid_t x = H5Dcreate(f, "vector", H5T_NATIVE_DOUBLE, s,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(x, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v);

    H5Dclose(x);
    H5Sclose(s);
    H5Fclose(f);
}

void load_vec_hdf5(const char *filename, double *v)
{
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t dataset_id = H5Dopen(file_id, "vector", H5P_DEFAULT);

    H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
            H5P_DEFAULT, v);

    H5Dclose(dataset_id);
    H5Fclose(file_id);
}


void build_random_sym_matrix(double* buffer,int n,double a ,double b )
{
   printf("Inside inocuous function\n");
    for (int i=0;i<n;i++)
    {
        buffer[matind(i,i,n)]=a+(b-a)*rand()*1./RAND_MAX;
        for (int j=0;j<i;j++)
        {
            buffer[matind(i,j,n)]=a+(b-a)*rand()*1./RAND_MAX;
            buffer[matind(j,i,n)]=buffer[matind(i,j,n)];
        }
    }
    return;
}

void print_matrix(double* buffer,int n)
{
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            printf("% 2.4f  ",buffer[matind(i,j,n)]);
        }
        printf("\n");
    }
    return;
}



void file_print_matrix(FILE* file, double* buffer,int n)
{
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            fprintf(file,"% 2.8f , ",buffer[matind(i,j,n)]);
        }
        fprintf(file,"\n");
    }
    return;
}

// void file_load_matrix()



void print_vector(double* buffer,int n)
{
    for (int i=0;i<n;i++)
    {
        printf("% 2.7f  ",buffer[i]); 
    }
    printf("\n");
    return;
}

void file_print_vector(FILE* file, double* buffer,int n)
{
    for (int i=0;i<n;i++)
    {
        fprintf(file,"% 2.7f",buffer[i]); 
        if (i<n-1)
        {
            fprintf(file," , "); 
        }
    }
    fprintf(file,"\n");
    return;
}

void householder(double* A,double*Z, double *d,double* e,int n)
{
    double* W =malloc(n*sizeof(double));
    double* V =malloc(n*sizeof(double));
    double* Q =malloc(n*sizeof(double));
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            Z[matind(i,j,n)]= (j==i) ? (1) : (0);
        }
    }
    for (int k=0;k<n-2;k++)
    {
        double s=formW(W,A,k,n);
        double c=formV(V,A,W,n);
        formQ(Q,V,W,c,n);
        formA(A,Q,W,n);
        formZ(Z,W,n);
        d[k]=A[matind(k,k,n)];
        e[k]=A[matind(k,k+1,n)];
    }
    d[n-2]=A[matind(n-2,n-2,n)];
    d[n-1]=A[matind(n-1,n-1,n)];
    e[n-2]=A[matind(n-2,n-1,n)];
    e[n-1]=0;
    free(W);
    free(V);
    free(Q);
}


void formZ(double* Z, double* W, int n)
{
    for (int i=0;i<n;i++)
    {
        double zeta=0;
        for (int k=0; k<n;k++)
        {
            zeta += Z[matind(i,k,n)]*W[k];
        }
        for (int j=0;j<n;j++)
        {
            Z[matind(i,j,n)] -= 2* zeta*W[j];
        }
        
    }
}

double formW(double* W, double* A, int k, int n)
{
    double s_2 = 0;
    for (int j=k+1;j<n;j++)
    {
        s_2+=A[matind(j,k,n)]*A[matind(j,k,n)];
    }
    double sign = (A[matind(k+1,k,n)] >= 0.) ? 1. : -1.;
    double s= sign*sqrt(s_2);
    double R=sqrt(2*s*(A[matind(k+1,k,n)]+s));
    W[k]=0;
    W[k+1]=(A[matind(k+1,k,n)]+s)/R;
    for (int j=k+2;j<n;j++)
    {
        W[j]=A[matind(j,k,n)]/R;
    }
    return s;
}

double formV(double* V, double* A, double* W, int n)
{
    double c=0;
    for (int i=0;i<n;i++)
    {
        V[i]=0;
        for (int j=0;j<n;j++)
        {
            V[i]+=A[matind(i,j,n)]*W[j];
        }
        c+= V[i]*W[i];
    }
    return c;
}

void formQ(double* Q, double* V, double* W, double c, int n)
{
    for (int j=0;j<n;j++)
    {
        Q[j]=V[j]-c*W[j];
    }
}

void formA(double* A, double* Q, double* W, int n)
{
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            A[matind(i,j,n)]-= 2*W[i]*Q[j] + 2*W[j]*Q[i];
        }
    }
}


void calc_block_start_end(int*block_start,int* block_end, double tol ,double* e,int n)
{
    *block_start=-1;
    for (int i=0; i<n-1;i++)
    {
        if (fabs(e[i])>=tol) 
        {
            *block_start=i;
            break;
        }
    }
    for (int i=*block_start+1; i<n;i++)
    {
        if (fabs(e[i])<=tol) 
        {
            *block_end=i;
            break;
        }
    }
}

void QL_iteration(double* Z, double*d,double* e, int n, double tol,int max_iter)
{
    int ql_iter=0;
    int block_start=0;
    int block_end=n-1;
    while (block_start >= 0 )
    {
        QL_iterate(Z,block_start,block_end,d,e,n);
        calc_block_start_end(&block_start,&block_end,tol,e,n);
        ql_iter++;
    }
    // printf("QL iterations = %d\n",ql_iter);
}


void QL_iterate(double*Z,int block_start, int block_end, double* d, double*e, int n)
{
    double shift;
    double tr=d[block_start]+d[block_start+1];
    double det=d[block_start]*d[block_start+1]-e[block_start]*e[block_start];
    // If block is 2x2, immediately calculate the eigenvalues
    if (block_end==block_start+1)
    {
        // printf("2x2 block \n");
        double mu = d[block_start]-d[block_start+1];
        double N=sqrt(mu*mu+4*e[block_start]*e[block_start]);
        double c=sqrt((1+mu/N)/2);
        double s=sqrt((1-mu/N)/2) * (e[block_start]>0 ? -1 : 1) ;
        e[block_start]=0;
        d[block_start]  =(tr+sqrt(tr*tr-4*det))/2;
        d[block_start+1]=(tr-sqrt(tr*tr-4*det))/2;
        for (int j=0;j<n;j++)
        {
            double F=Z[matind(j,block_start,n)];
            Z[matind(j,block_start  ,n)]=c*F-s*Z[matind(j,block_start+1,n)];
            Z[matind(j,block_start+1,n)]=s*F+c*Z[matind(j,block_start+1,n)];
        }
        return;
    }
    
    if (d[block_start]>tr/2)
    {
        shift=(tr+sqrt(tr*tr-4*det))/2;
    }
    else
    {
        shift=(tr-sqrt(tr*tr-4*det))/2;
    }

    
    // Otherwise, continue and perform QL iteration on the first block
    // Replace n-1 by block_end and 0 by block_start
    shift=0;
    double p= d[block_end] - shift ;
    double c=1;
    double s=0;
    double crit=0;
    
    for (int i=block_end-1;i>=block_start;i--)
    {
        double r=sqrt(p*p+e[i]*e[i]);
        double g=c*e[i];
        double h=c*p;
        e[i+1]=s*r;
        c=p/r;
        s=e[i]/r;
        p=c*(d[i]-shift)-s*g;
        d[i+1]=h+s*(c*g+s*(d[i]-shift))+shift;
        for (int j=0;j<n;j++)
        {
            double F=Z[matind(j,i,n)];
            Z[matind(j,i  ,n)]=c*F-s*Z[matind(j,i+1,n)];
            Z[matind(j,i+1,n)]=s*F+c*Z[matind(j,i+1,n)];
        }
    }
    e[block_start]=s*p;
    d[block_start]=c*p+shift;
}