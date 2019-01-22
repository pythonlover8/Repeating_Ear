#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>
#include "model_search.h"
#include "matrix.h"
#include "nrutil.h"
#include "weight.h"
/* conduct bootstrap analysis*/
/*input: the file which stores the xcorr results ,repeating time, evdp,  
output: the file which stores the Empirical probability distributions of the model parameters. */
int main(int argc, char **argv){
	if(argc!=4){
		nrerror("Please input 3 arguments: the file which stores the xcorr results ,repeating time, evdp");
	}
	int repeat_time = atof(argv[2]),i,ir,j,iter=1,Niter=2;
	int **index,count,M=4;
	float evdp = atof(argv[3]),a_radius=6371.1;
	float **result,**result_resampled,**result_resampled1;
	float **W,**A,**B,**D,**F,**invA,E=0.0,**G,**GT;
	float *mest,*dobs,*C,*e;
	char filename[160],file[160];
	strcpy(file,argv[1]);
	FILE *fp = fopen(file,"r"), *fp1 = fopen("model.dat","w");
	if(fp==NULL||fp1==NULL){
		nrerror("Unable to open the file");
	}

	source sour;
	strcpy(filename,"/Users/a123/Desktop/Ear_relo/iasp91.dat");
	sour = model_search(evdp,a_radius,filename);
	result = matrix(1,MAX_Size,1,5);
	for(i=1;!feof(fp);i++){
		fscanf(fp,"%f %f %f %f %f", &result[i][1], &result[i][2], &result[i][3],
			 &result[i][4], &result[i][5]);
	}
	fclose(fp);

	count = i - 2;
	mest = vector(1,M);
	dobs = vector(1,count);
	C = vector(1,count);
	e = vector(1,count);
	A = matrix(1,M,1,count);
	B = matrix(1,count,1,count);
	D = matrix(1,M,1,M);
	F = matrix(1,count,1,M);
	invA = matrix(1,M,1,M);
	G = matrix(1,count,1,M);
	GT = matrix(1,M,1,count);
	result_resampled = matrix(1,5,1,count);
	result_resampled1 = matrix(1,count,1,5);
	for(ir=1;ir<=repeat_time;ir++){
		if(ir==1){
			for(i=1;i<=count;i++){
				G[i][1] = -1/sour.vp*sin(result[i][5]*pi/180)*cos(result[i][2]*pi/180);
				G[i][2] = -1/sour.vp*sin(result[i][5]*pi/180)*sin(result[i][2]*pi/180);
				G[i][3] = 1/sour.vp*cos(result[i][5]*pi/180);
				G[i][4] = 1;
				dobs[i]=result[i][3];
			}
		}	
		else{
		
		index = unirand(1,count,count,1);
		result_resampled[1] = get_vector_mat(result,index,count,1);
		result_resampled[2] = get_vector_mat(result,index,count,2);  
		result_resampled[3] = get_vector_mat(result,index,count,3);
		result_resampled[4] = get_vector_mat(result,index,count,4);
		result_resampled[5] = get_vector_mat(result,index,count,5);
		transpose(result_resampled,result_resampled1,5,count);

		for(i=1;i<=count;i++){
				G[i][1] = -1/sour.vp*sin(result_resampled1[i][5]*pi/180)*cos(result_resampled1[i][2]*pi/180);
				G[i][2] = -1/sour.vp*sin(result_resampled1[i][5]*pi/180)*sin(result_resampled1[i][2]*pi/180);
				G[i][3] = 1/sour.vp*cos(result_resampled1[i][5]*pi/180);
				G[i][4] = 1;
				dobs[i]=result_resampled1[i][3];
			}
		}
		/* this is the block for downweighting the time residuals iteratively */
			
			for(iter=1;iter<=Niter;iter++){
				if(iter==1){
					W = eye(count);

				}
				else{
					weight(e,count,W);
				}
				transpose(G,GT,count,M);
				matmul(GT,M,count,W,count,A);
				matmul(A,M,count,G,M,D);
				inv(D,M);
				matmul(D,M,M,GT,count,B);
				matmul(B,M,count,W,count,B);
				matmul2(B,M,count,dobs,mest);

				matmul2(G,count,M,mest,C);
				for(i=1;i<=count;i++){
					e[i]=dobs[i]-C[i];
				}
			}
			matmul3(e,W,count,count,e);
			E=dot_product1(e,e,count);
			E = sqrt(E/(count-M));
		fprintf(fp1,"%7.4f %7.4f %7.4f %7.4f %8.6f\n",mest[1],mest[2],
			mest[3],mest[4],E);
	}
	fclose(fp1);
	free_vector(e,1,count);
	free_vector(mest,1,M);
	free_vector(dobs,1,count);
	free_vector(C,1,count);
	free_imatrix(index,1,count,1,1);
	free_matrix(A,1,count,1,count);
	free_matrix(W,1,count,1,count);
	free_matrix(B,1,count,1,count);
	free_matrix(invA,1,M,1,M);
	free_matrix(G,1,count,1,M);
	free_matrix(GT,1,M,1,count);
	free_matrix(result,1,MAX_Size,1,5);
	free_matrix(result_resampled,1,5,1,count);
	free_matrix(result_resampled1,1,count,1,5);
	return 0;
}