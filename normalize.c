#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>
#include "model_search.h"
#include "sac.h"
#include "matrix.h"
#include "nrutil.h"

/*input: the list of sacfiles, scaling factor, and time-window */
int main(int argc,char **argv){
	int i,ii=1,j,jj,nn1,nn2;
	int **index=NULL;
	char sacfile[160],saclist[160];
	float scaling_factor,t1,t2;
	float *trace=NULL,*data=NULL,sign;
	max max1;
	min min1;
	SACHEAD hd;
	strcpy(saclist,argv[1]);
	scaling_factor=atof(argv[2]);
	t1=atof(argv[3]);
	t2=atof(argv[4]);
	FILE *f1;
	if((f1=fopen(saclist,"r"))==NULL){
		nrerror("unable to open the file");
	}
	for(i=1;!feof(f1);i++){
		sign=1.0;
		ii=1;
		fscanf(f1,"%s",sacfile);
		if((trace=read_sac(sacfile,&hd))==NULL){
			nrerror("unable to open the sacfile");
		}		
		nn1=floor((t1-hd.b)/hd.delta);
		nn2=ceil((t2-hd.b)/hd.delta);
		data = vector(1,nn2-nn1+1);
		index = imatrix(1,nn2-nn1+1,1,1);
		for(j=nn1;j<=nn2;j++){
			index[ii][1] = j;
			ii++;
		}
		data=get_vector_vec(trace,index,nn2-nn1+1);
		max1=max_vector(data,nn2-nn1+1);
		min1=min_vector(data,nn2-nn1+1);
	/*	if(fabs(max1.max_value)<fabs(min1.min_value)){
			sign=-1.0;
		}	*/
		for(jj=1;jj<=hd.npts;jj++){
			trace[jj]/=(max1.max_value-min1.min_value);
			trace[jj]=trace[jj]*sign*scaling_factor;
		}
		write_sac(sacfile,hd,trace);
		free_vector(data,1,nn2-nn1+1);
		free_imatrix(index,1,nn2-nn1+1,1,1);
	    free(trace);
	}
	fclose(f1);
	return 0;
}
