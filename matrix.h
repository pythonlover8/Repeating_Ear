#ifndef _MATRIX_H_
#define _MATRIX_H_
#include "nrutil.h"
#include "model_search.h"
#include "matrix.h"
void transpose(float **A, float **B, int N, int M);
float dot_product1(float *xarray, float *yarray, int count);
void matmul(float **A,int N,int M,float **B,int MB,float **C);
void matmul2(float **A,int N,int M,float *B,float *C);
void matmul3(float *A,float **B, int N, int M,float*C);
float **eye(int N);
void ludcmp(float **a, int n, int *indx, float *d);
void lubksb(float **a, int n, int *indx, float b[]);
void gaussj(float **a,int n,float **b, int m);
void inv(float **A,int N);
void print_mat(float **a,int n,int M);
void print_vec(float *a,int n);
#endif /*_MATRIX_H_*/	