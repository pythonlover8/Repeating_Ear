#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "matrix.h"
#include "model_search.h"
#include "nrutil.h"
#include "weight.h"

/*set up the weighting matrix, using a biweight function proposed by Mosteller and Tukey, 1977.
input: a vector of time residuals, size of the vector;
output: a weight matrix*/
void weight(float *e,int N,float **W){
	int i;
	float *dr,drmad,emedian;
	float *dt_wt;
	dr=vector(1,N);
	dt_wt = vector(1,N);
	emedian=median(e,N);
	for(i=1;i<=N;i++){
		dr[i] = fabs(e[i]-emedian);
	}
	drmad = median(dr,N);
	drmad = drmad*alpha/sigma_mad;
	for(i=1;i<=N;i++){
		dt_wt[i] = 1 - (e[i]/drmad)*(e[i]/drmad);
		if(dt_wt[i]<0){
			W[i][i]=0;
		}
		else{
			W[i][i]=(dt_wt[i])*(dt_wt[i]);
		}
	}
	free_vector(dt_wt,1,N);
	free_vector(dr,1,N);

}