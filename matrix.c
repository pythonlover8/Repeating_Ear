#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "nrutil.h"
#include "model_search.h"
#include "matrix.h"
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define TINY 1.0e-20;
/* ludcmp.c,lubksb.c, and gaussj.c are scripts written by Press, W.H., Teukolsky, S.A., Vetterling, W.T., & Flannery, B.P. (1997)*/
void transpose(float **A, float **B, int N, int M){
	int i;
	for (i=1;i<=M;i++){
		for(int j=1;j<=N;j++){
			B[i][j] = A[j][i];
		}
	}
}
float dot_product1(float *xarray, float *yarray, int count){
	int i;
	float sum=0;
	for(i=1;i<=count;i++){
		sum += xarray[i]*yarray[i];
	}
	return sum;
}
/* for the case: A(n,m)*B(m,mb)=C(n,mb);*/
void matmul(float **A,int N,int M,float **B,int MB,float **C){
	int i,j,k;
	float sum;
	for(i=1;i<=N;i++){
		for(j=1;j<=MB;j++){
			sum=0.0;
			for(k=1;k<=M;k++)
			{		
				sum += A[i][k]*B[k][j];
		}
		C[i][j]=sum;
	}
}
}
/* for the A(n,m)*B(m) = C(n);*/
void matmul2(float **A,int N,int M,float *B,float *C){
	int i,j;
	float sum;
	for(i=1;i<=N;i++){
		sum=0.0;
		for(j=1;j<=M;j++){
			sum += A[i][j]*B[j];
		}
		C[i]=sum;

	}
}
/* for the A(n)*B(n,m) = C(m);*/
void matmul3(float *A,float **B, int N, int M,float*C){
	int i,j;
	float sum;
	for(i=1;i<=M;i++){
		sum=0.0;
		for(j=1;j<=N;j++){
			sum += A[j]*B[j][i];
		}
		C[i]=sum;
	}
}
float **eye(int N){
	float **A;
	int i;
	A = matrix(1,N,1,N);
	for(i=1;i<=N;i++){
		A[i][i] = 1.0;
	}
	return A;
}
void ludcmp(float **a, int n, int *indx, float *d)
{
	int i,imax,j,k;
	float big,dum,sum,temp;
	float *vv;

	vv=vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_vector(vv,1,n);
}


void lubksb(float **a, int n, int *indx, float b[])
{
	int i,ii=0,ip,j;
	float sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
void gaussj(float **a,int n,float **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	float big,dum,pivinv,temp;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l]);
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l]);
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}
void inv(float **A,int N){
	float **B;
	int i;
	B = eye(N);
	gaussj(A,N,B,N);

	
	free_matrix(B,1,N,1,N);
}

void print_mat(float **a,int n,int M){
	int i,j;
	for(i=1;i <= n;i++){
		for(j=1;j <= M;j++)
			printf(" %.2f ",a[i][j]);
		printf("\n");
		}
}

void print_vec(float *a,int n){
	int i,j;
	for(i=1;i <= n;i++){
		
	    printf(" %7.4f \t",a[i]);
		printf("\n");
		}
}
