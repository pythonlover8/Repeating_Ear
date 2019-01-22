#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "sac.h"
#include "model_search.h"
#include "nrutil.h"

typedef struct coeff
{
	float max_c;
	float time_lag;
}coeff;

float dot_product1(float *xarray, float *yarray, int count){
	int i;
	float sum=0;
	for(i=1;i<=count;i++){
		sum += xarray[i]*yarray[i];
	}
	return sum;
}


coeff xcorr(float *trace1, float *trace2,float dt,int count){
	int i,ii,m;
	float *sum_c,*psi,sum_x2,sum_y2;
	float aa,time_lag;
	coeff c;
	max psi_max;
	sum_x2 = 0;
	sum_y2 = 0;
	sum_c = vector(1,2*count-1);
	psi = vector(1,2*count-1);
	sum_x2 = dot_product1(trace1,trace1,count);
	sum_y2 = dot_product1(trace2,trace2,count);
	aa = sqrt(sum_y2)/sqrt(sum_x2);
	for(i=1;i<=2*count-1;i++){
		sum_c[i] = 0;
	}
	ii = 1;
	for(i=1-count;i<count;i++){
		if(i<=0){
			for(m=1;m<=count+i;m++){
				sum_c[ii] += trace1[m]*trace2[m-i];
			}
		}
		else
		{
			for(m=1;m<=count-i;m++){
			    sum_c[ii] += trace1[m+i]*trace2[m];
			}
		}
		psi[ii] = sum_c[ii]*aa/sum_y2;
		ii++;
	}
	psi_max = max_vector(psi,2*count-1);
	c.max_c = psi_max.max_value;
	c.time_lag = (psi_max.max_loc-count+1)*dt;
	free_vector(psi,1,2*count-1);
	free_vector(sum_c,1,2*count-1);
	return c;
}



/* first input the file where stores the names of the sacfiles to be processes
then the template file, finally the timewindow t1,t2. */
/* do the resampling first*/
int main(int argc, char **argv){

	int i,ii,count,n1,nn1,n2,nn2;
	int **index;
	char filename[160],tempfile[160];
	float *trace1,*trace2,*data,*data1,dt,dtemp;
	float t1,t2;
	coeff c;
	if(argc<5) nrerror("Too few arguments supplied!");
	if(argc>5) nrerror("Too many arguments supplied!");
	FILE *f,*f1=fopen("unshiftedd.dat","w");
	SACHEAD hd1,hd2;

	if((f=fopen(argv[1],"r"))==NULL||f1==NULL){
		nrerror("Unable to open the file");
	}
	t1 = atof(argv[3]);
	t2 = atof(argv[4]);
	strcpy(tempfile,argv[2]);

	if((trace2=read_sac(tempfile,&hd2))==NULL){
		nrerror("Unable to open the sacfile");
	}
	n2 = floor((t1-hd2.b)/hd2.delta);
	nn2 = ceil((t2-hd2.b)/hd2.delta);
	index = imatrix(1,nn2-n2+1,1,1);
	ii = 1;
	for(i=n2;i<=nn2;i++){
		index[ii][1] = i;
		ii++;
	}
	data = vector(1,nn2-n2+1);
	data = get_vector_vec(trace2,index,nn2-n2+1);
	while(!feof(f)){
	   fscanf(f,"%s",filename); 
	   if((trace1 = read_sac(filename,&hd1))==NULL){
	   	   nrerror("Unable to open the sacfile");
	   }
	   if(hd1.delta!=hd2.delta){
	   	 printf("Time_intervals of two files are not equal\n");
	   	 printf("Not processing %s ...\n",filename);
	   	 continue;
	   }
	   n1 = floor((t1-hd1.b)/hd2.delta);
	   nn1 = ceil((t2-hd1.b)/hd2.delta);
	   data1 = vector(1,nn1-n1+1);
	   ii = 1;
	   for(i=n1;i<=nn1;i++){
		index[ii][1] = i;
		ii++;
	   }
	   data1 = get_vector_vec(trace1,index,nn1-n1+1);
	   c = xcorr(data,data1,hd2.delta,nn1-n1+1);
	   if(c.max_c<0.60){
	   	printf("max coeff = %f, which is below the lowest level\n do not shift the seismogram...\n",c.max_c);
	   	fprintf(f1,"%s\n",filename);
	   	continue;
	   }
	   hd1.b = hd1.b+c.time_lag;
	   printf("%f %f %f %f\n",hd1.b-c.time_lag,c.time_lag,hd1.b,c.max_c);
	   write_sac(filename, hd1, trace1);
	   free_vector(data1,1,nn1-n1+1);
	   free(trace1);
}
	
	free_vector(data,1,nn2-n2+1);
	fclose(f);
	fclose(f1);
	free_imatrix(index,1,nn2-n2+1,1,1);
	
	free(trace2);
	return 0;
}
