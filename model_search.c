#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "nrutil.h"
#include"model_search.h"

source model_search(float evdp, float a_radius, char *file){
	int i;
	float *radius=NULL,*alpha_s=NULL,*beta_s=NULL;
	float *z_s=NULL;
	source sour;
	FILE *f;
	radius = vector(1,MAX_Array);
	alpha_s = vector(1,MAX_Array);
	beta_s = vector(1,MAX_Array);
	z_s = vector(1,MAX_Array);
    f = fopen(file,"r");
	if (f==NULL)
	{
		nrerror("unable to open the file");
	}
    sour.vp = 0;
    sour.ev_radius = 0;
	for(i=1;!feof(f);i++){
		fscanf(f,"%f %f %f",&radius[i],&alpha_s[i],&beta_s[i]);
		z_s[i] = a_radius - radius[i];
		if(i==1) continue;
		if(evdp>=z_s[i-1]&&evdp<=z_s[i]){
			sour.vp = alpha_s[i-1]+(evdp-z_s[i-1])/(z_s[i]-z_s[i-1])
				*(alpha_s[i]-alpha_s[i-1]);
			sour.ev_radius = a_radius-evdp;
			printf("vp = %f\n",sour.vp);
			printf("ev_radius = %f\n",sour.ev_radius);
			break;
		}
	} 
	fclose(f);	
	return sour;
	free_vector(radius,1,MAX_Array);
	free_vector(alpha_s,1,MAX_Array);
	free_vector(beta_s,1,MAX_Array);
	free_vector(z_s,1,MAX_Array);
}

float* linespace(float fl, float fr, int n){
	int i;
    float *u,interval;
    u = vector(1,n);
	interval = (fr-fl)/(n-1);
	for(i=1;i<=n-1;i++){
		u[i] = fl + (i-1) * interval;
	}
	u[n] = fr;
	return u;
}

float* range(float fl, float fr, float interval){
	int i,n;
	float *u;
	n = ceil((fr-fl)/interval);
	u = vector(1,n);
	for(i=1;i<=n;i++){
		u[i] = fl + (i-1)*interval;
	}
	return u;
}
int** unirand(int min_num, int max_num, int N, int M){
	int i,j;
	int **rand_num;
	rand_num = imatrix(1,N,1,M);
	for(i=1;i<=N;i++){
		for(j=1;j<=M;j++){
			rand_num[i][j] = rand() % (max_num+1-min_num) + min_num;
 		}
	}
	return rand_num;
}
float* getcol(float **mat,int col){
	int nrows = sizeof(mat)/sizeof(mat[0]);
	int ncols = sizeof(mat[0])/sizeof(mat[0][0]);
	int sz = nrows*ncols,n=0;
	float *out;
	out = vector(1,nrows);
	for(int i=1;i<=sz;i+=ncols){
		out[n] = mat[i][col];
		n++;
	}
	return out;
} 
float* get_vector_mat(float **mat, int **index, int count, int col){
	float *out;
	int i,nind;
	out = vector(1,count);
	for(i=1;i<=count;i++){
		nind = index[i][1];
		out[i] = mat[nind][col];
	}
	return out;
} 

float* get_vector_vec(float *vec, int **index, int count){
	float *out;
	int i,nind;
	out = vector(1,count);
	for(i=1;i<=count;i++){
		nind = index[i][1];
		out[i] = vec[nind];
	}
	return out;
} 
max max_vector(float *v, int count){
	int i;
	max m;
	m.max_value = -100000;
	m.max_loc = 0;
	for(i=1;i<=count;i++){
		if(v[i]>m.max_value){
			m.max_value = v[i];
			m.max_loc = i;
		}
	}
	return m;
}
min min_vector(float *v, int count){
	int i;
	min m;
	m.min_value = 100000;
	m.min_loc = 0;
	for(i=1;i<=count;i++){
		if(v[i]<m.min_value){
			m.min_value = v[i];
			m.min_loc = i;
		}
	}
	return m;
}
void swap_float(float xp,float yp){
	float temp=xp;
	xp=yp;
	yp=temp;
}
void swap_int(int xp,int yp){
	int temp=xp;
	xp=yp;
	yp=temp;
}
/*sort the vector float elements in ascending order
using the selection sort algorithm*/
void sort_float(float *v,int N){
	int i,j,smallNdx;
	for(i=1;i<=N;i++){
		smallNdx=i;
		for(j=i+1;j<=N;j++){
			if(v[j]<v[smallNdx]){
				smallNdx=j;
			}
		}
		if(smallNdx!=i){
			swap_float(v[j],v[smallNdx]);
		}
	}
}
void sort_int(int *v,int N){
	int i,j,smallNdx;
	for(i=1;i<=N;i++){
		smallNdx=i;
		for(j=i+1;j<=N;j++){
			if(v[j]<v[smallNdx]){
				smallNdx=j;
			}
		}
		if(smallNdx!=i){
			swap_int(v[j],v[smallNdx]);
		}
	}
}
/* find median value of a vector*/
float median(float x[], int n) {
    float temp;
    int i, j;
    // the following two loops sort the array x in ascending order
    for(i=1; i<=n-1; i++) {
        for(j=i+1; j<=n; j++) {
            if(x[j] < x[i]) {
                // swap elements
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }
    	
    }

    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[n/2] + x[n/2 + 1]) / 2.0);
    } else {
        // else return the element in the middle
        return x[n/2+1];
    }
}
