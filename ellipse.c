#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>
#include<time.h>
#include "model_search.h"
#include "nrutil.h"
#include "matrix.h"

bool if_same_qua(int quadrant,float coordinate1,float coordinate2){
	if(quadrant==1){
		if(coordinate1>0&&coordinate2>0){
			return true;
		}
		else
		{
			return false;
		}
	}
	else if(quadrant==2){
		if(coordinate1>0&&coordinate2<0){
			return true;
		}
		else{
			return false;
		}
	}
	else if(quadrant==3){
		if(coordinate1<0&&coordinate2<0){
			return true;
		}
		else{
			return false;
		}
	}
	else{
		if(coordinate1<0&&coordinate2>0){
			return true;
		}
		else{
			return false;
		}
	}
}
/* input the file that stores the results
output the error ellispe*/
int main(int argc, char **argv){
	int i,j,count,*quadrant,sign=1;
	if(argc<1){
		nrerror("Please input: file that stores the results,\n");
	}
	char filename1[160];
	strcpy(filename1,argv[1]);
	FILE *fp1=fopen(filename1,"r"),*fp2=fopen("points.dat","w");
	FILE *fp3=fopen("ellipse.dat","w");
	float **result,*strike,**cov,*eigval,*e;
	float *hsep,mean_hsep,mean_depth,a,b;
	float chisquare_val=2.4477,angle;
	float max_eig,min_eig,*max_evc,*min_evc;
	result=matrix(1,5,1,MAX_Array);
	cov=matrix(1,2,1,2);
	strike=vector(1,2);
	eigval=vector(1,2);
	e=vector(1,2);
	max_evc=vector(1,2);
	min_evc=vector(1,2);
	for(i=1;!feof(fp1);i++){
		fscanf(fp1,"%f %f %f %f %f",
			&result[1][i],&result[2][i],&result[3][i],&result[4][i],&result[5][i]);
	}
	count = i-2;
	fclose(fp1);
	hsep=vector(1,count);
	quadrant=ivector(1,2);
	printf("Input: strikes of two possible fault planes\n");
	scanf("%f %f",&strike[1],&strike[2]);
	quadrant[1]=ceil(strike[1]/90.0);
	quadrant[2]=ceil(strike[2]/90.0);
	for(i=1;i<=count;i++){
		if(if_same_qua(quadrant[1],result[1][i],result[2][i])||
			if_same_qua(quadrant[2],result[1][i],result[2][i]))
		{
			sign=1;
		}
		else{
			sign=-1;
		}
		hsep[i]=sign*sqrt(result[1][i]*result[1][i]+
			result[2][i]*result[2][i]);	
		fprintf(fp2,"%f %f\n",hsep[i],result[3][i]);
	}
	fclose(fp2);
	/* fit the error ellipse*/
	/* compute the covariance matrix */
	mean_hsep=mean(hsep,count);
	mean_depth=mean(result[3],count);
	printf("mean_hsep=%f\n",mean_hsep);
	printf("mean_depth=%f\n",mean_depth);
	for(i=1;i<=count;i++){
		/*fill in the first entry*/
		cov[1][1]+=(hsep[i]-mean_hsep)*(hsep[i]-mean_hsep);
		/*fill in the second entry*/
		cov[1][2]+=(hsep[i]-mean_hsep)*(result[3][i]-mean_depth);
		/*fill in the third entry*/
		cov[2][1]+=(result[3][i]-mean_depth)*(hsep[i]-mean_hsep);
		/*fill in the fourth entry*/
		cov[2][2]+=(result[3][i]-mean_depth)*(result[3][i]-mean_depth);
	}
	cov[1][1]/=(count-2);
	cov[1][2]/=(count-2);
	cov[2][1]/=(count-2);
	cov[2][2]/=(count-2);

	print_mat(cov,2,2);
	tred2(cov,2,eigval,e);
	tqli(eigval,e,2,cov);
	print_mat(cov,2,2);
	max_evc=getcol(cov,1,2,2);
	print_vec(max_evc,2);
	max_eig=eigval[1];
	min_evc=getcol(cov,2,2,2);
	min_eig=eigval[2];
	angle=atan(max_evc[2]/max_evc[1]);
	angle=angle*180/pi;
	if(angle<0){
		angle=angle+360;
	}
	a=chisquare_val*sqrt(max_eig);
	b=chisquare_val*sqrt(min_eig);
	printf("major axis=%f, minor axis=%f, angle=%f\n",a,b,angle);
	printf("center_x=%f,center_y=%f\n",mean_hsep,mean_depth);
	fprintf(fp3,"%f %f %f %f %f\n",a,b,angle,mean_hsep,mean_depth);
	fclose(fp3);
	free_vector(strike,1,2);
	free_vector(eigval,1,2);
	free_vector(e,1,2);
	free_vector(max_evc,1,2);
	free_vector(min_evc,1,2);
	free_vector(hsep,1,count);
	free_ivector(quadrant,1,2);
	free_matrix(result,1,5,1,MAX_Array);
	free_matrix(cov,1,2,1,2);
	return 0;
}
