#ifndef EIGEN_SERIAL_H
#define EIGEN_SERIAL_H

#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>

#define matind(i,j,n) (((i)*(n))+(j))

void save_vec_hdf5(const char *file, double *v, int len);
void load_vec_hdf5(const char *filename, double *v);

void build_random_sym_matrix(double* ,int,double,double);

void print_matrix(double*,int);
void print_vector(double* buffer,int n);
void file_print_matrix(FILE* file,double* buffer,int n);
void file_print_vector(FILE* file,double* buffer,int n);



void householder(double* A,double* Z,double* d,double* e,int n);


double formW(double* W, double* A, int k, int n);
double formV(double* V, double* A, double* W, int n);
void formQ(double* Q, double* V, double* W, double c, int n);
void formA(double* A, double* Q, double* W, int n);
void formZ(double* Z, double* W, int n);

void calc_block_start_end(int*block_start,int* block_end, double tol ,double* e,int n);
void QL_iterate(double*Z,int block_start, int block_end, double* d, double*e, int n);
void QL_iteration(double*Z, double*d, double* e,int n , double tol,int max_iter);
#endif