#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h>
#include "model_search.h"
#include "nrutil.h"
#include "matrix.h"

/* input file, the no. of bins */
int main(int argc, char **argv){
	int i,j,count, Nbins,M=4, *fre,sum=0;
	if(argc<4){
		nrerror("Please input: file, the no. of bins, mode \n");
	}
	Nbins = atof(argv[2]);
	char filename[160];
	float **result,normalize=atof(argv[3]);
	float *m1bins,*m2bins,*m3bins,*m4bins,*pm,*Pm;
	float Dm1bins,Dm2bins,Dm3bins,Dm4bins,m1low,m1high;
	float m2low,m2high,m3low,m3high,m4low,m4high,*mest;
	min min1,min2,min3,min4;
	max max1,max2,max3,max4,maxP1,maxP2,maxP3,maxP4;
	bool flag,flag1;
	strcpy(filename,argv[1]);
	FILE *fp=fopen(filename,"r"),*fp1=fopen("prob1.dat","w");
	FILE *fp2=fopen("prob2.dat","w"),*fp3=fopen("prob3.dat","w");
	FILE *fp4=fopen("prob4.dat","w");
	FILE *fp5=fopen("hist_interval.dat","w"),*fp6=fopen("maxmin.dat","w");
	FILE *fp7=fopen("range.dat","w"),*fp8=fopen("mest.dat","w");
	FILE *fp9=fopen("interval.dat","w");
	if(fp==NULL||fp1==NULL||fp2==NULL||fp3==NULL||fp4==NULL
		||fp5==NULL||fp6==NULL||fp7==NULL||fp8==NULL){
		nrerror("Unable to open the file");
	}
	result = matrix(1,5,1,MAX_Array);
	for(i=1;!feof(fp);i++){
		fscanf(fp,"%f %f %f %f %f",
			&result[1][i],&result[2][i],&result[3][i],&result[4][i],&result[5][i]);
	}
    count = i-2;
    fclose(fp);
    mest=vector(1,M);
	min1 = min_vector(result[1],count);
	max1 = max_vector(result[1],count);
	Dm1bins=(max1.max_value-min1.min_value)/(Nbins-1);
	m1bins = linespace(min1.min_value,max1.max_value,Nbins);
	min2 = min_vector(result[2],count);
	max2 = max_vector(result[2],count);
	Dm2bins=(max2.max_value-min2.min_value)/(Nbins-1);
	m2bins = linespace(min2.min_value,max2.max_value,Nbins);
	min3 = min_vector(result[3],count);
	max3 = max_vector(result[3],count);
	Dm3bins=(max3.max_value-min3.min_value)/(Nbins-1);
	m3bins = linespace(min3.min_value,max3.max_value,Nbins);
	min4 = min_vector(result[4],count);
	max4 = max_vector(result[4],count);
	Dm4bins=(max4.max_value-min4.min_value)/(Nbins-1);
	m4bins = linespace(min4.min_value,max4.max_value,Nbins);
	fre = ivector(1,Nbins);
	pm = vector(1,Nbins);
	Pm = vector(1,Nbins);
	for(i=1;i<Nbins;i++){
		fre[i] = 0;
		for(j=1;j<=count;j++){
			if(result[1][j]<m1bins[i+1]&&result[1][j]>=m1bins[i]){
				fre[i] += 1;
			}
		}
	}	

	for(i=1;i<Nbins;i++){
		sum += fre[i];
	}
	flag=true;
	flag1=true;
	if(normalize){
	for(i=1;i<Nbins;i++){
		if(i==1){
		pm[i] = (float) fre[i]/sum;
		Pm[i] =  pm[i];
		}
		else{
			pm[i] = (float) fre[i]/sum;
			Pm[i] = Pm[i-1] + pm[i];
		}
		if(Pm[i]>=0.05&&flag) {
			fprintf(fp7, "%f %f %f\n",
			(m1bins[i+1]+m1bins[i])/2.0,pm[i],Pm[i]); 
			flag=false;
		}
		if(Pm[i]>=0.95&&flag1){
			fprintf(fp7, "%f %f %f\n",(m1bins[i+1]+m1bins[i])/2.0,pm[i],Pm[i]);
			flag1=false;
		}
		fprintf(fp1,"%7.4f %6.4f %6.4f\n", (m1bins[i+1]+m1bins[i])/2.0,pm[i],Pm[i]);
	}
}
	else{
		for(i=1;i<Nbins;i++){
		if(i==1){
		pm[i] = (float) fre[i]/sum;
		Pm[i] =  pm[i];
		}
		else{
			pm[i] = (float) fre[i]/sum;
			Pm[i] = Pm[i-1] + pm[i];
		}
		if(Pm[i]>=0.05&&flag) {
			fprintf(fp7, "%f %f %f\n",
			(m1bins[i+1]+m1bins[i])/2.0,pm[i],Pm[i]); 
			flag=false;
		}
		if(Pm[i]>=0.95&&flag1){
			fprintf(fp7, "%f %f %f\n",(m1bins[i+1]+m1bins[i])/2.0,pm[i],Pm[i]);
			flag1=false;
		}
		fprintf(fp1,"%7.4f %6.4f %6.4f\n", (m1bins[i+1]+m1bins[i])/2.0,pm[i]*sum,Pm[i]*sum);
	}
}
	fclose(fp1);
	maxP1=max_vector(pm,Nbins-1);
	if(normalize) {
		fprintf(fp8,"%f %f\n",(maxP1.max_loc-1)*Dm1bins+(m1bins[2]+m1bins[1])/2.0,
		maxP1.max_value);
		fprintf(fp9,"%d\n",(int)(maxP1.max_value)/4);
	}
	else {
		fprintf(fp8,"%f %f\n",(maxP1.max_loc-1)*Dm1bins+(m1bins[2]+m1bins[1])/2.0,
		maxP1.max_value*sum);
		fprintf(fp9,"%d\n",(int)(maxP1.max_value*sum)/4);
	}
	for(i=1;i<Nbins;i++){
		fre[i] = 0;
		for(j=1;j<=count;j++){
			if(result[2][j]<m2bins[i+1]&&result[2][j]>=m2bins[i]){
				fre[i] += 1;
			}
		}
	}	
	sum=0;
	for(i=1;i<Nbins;i++){
		sum += fre[i];
	}
	flag=true;
	flag1=true;
	if(normalize){
	for(i=1;i<Nbins;i++){
		if(i==1){
		pm[i] = (float) fre[i]/sum;

		Pm[i] =  pm[i];
		}
		else{
			pm[i] = (float) fre[i]/sum;
			Pm[i] = Pm[i-1] + pm[i];
		}
		if(Pm[i]>=0.05&&flag) {
			fprintf(fp7, "%f %f %f\n",
			(m2bins[i+1]+m2bins[i])/2.0,pm[i],Pm[i]); 
			flag=false;
		}
		if(Pm[i]>=0.95&&flag1){
			fprintf(fp7, "%f %f %f\n",(m2bins[i+1]+m2bins[i])/2.0,pm[i],Pm[i]);
			flag1=false;
		}
		fprintf(fp2,"%7.4f %6.4f %6.4f\n", (m2bins[i+1]+m2bins[i])/2.0,pm[i],Pm[i]);
	}
}
	else{
		for(i=1;i<Nbins;i++){
		if(i==1){
		pm[i] = (float) fre[i]/sum;
		Pm[i] =  pm[i];
		}
		else{
			pm[i] = (float) fre[i]/sum;
			Pm[i] = Pm[i-1] + pm[i];
		}
		if(Pm[i]>=0.05&&flag) {
			fprintf(fp7, "%f %f %f\n",
			(m2bins[i+1]+m2bins[i])/2.0,pm[i],Pm[i]); 
			flag=false;
		}
		if(Pm[i]>=0.95&&flag1){
			fprintf(fp7, "%f %f %f\n",(m2bins[i+1]+m2bins[i])/2.0,pm[i],Pm[i]);
			flag1=false;
		}
		fprintf(fp2,"%7.4f %6.4f %6.4f\n", (m2bins[i+1]+m2bins[i])/2.0,pm[i]*sum,Pm[i]*sum);
	}
}
	fclose(fp2);
	maxP2=max_vector(pm,Nbins-1);
	if(normalize) {
		fprintf(fp8,"%f %f\n",(maxP2.max_loc-1)*Dm2bins+(m2bins[2]+m2bins[1])/2.0,
		maxP2.max_value);
		fprintf(fp9,"%d\n",(int)(maxP2.max_value)/4);
	}
	else {
		fprintf(fp8,"%f %f\n",(maxP2.max_loc-1)*Dm2bins+(m2bins[2]+m2bins[1])/2.0,
		maxP2.max_value*sum);
		fprintf(fp9,"%d\n",(int)(maxP2.max_value*sum)/4);
	}
		for(i=1;i<Nbins;i++){
		fre[i] = 0;
		for(j=1;j<=count;j++){
			if(result[3][j]<m3bins[i+1]&&result[3][j]>=m3bins[i]){
				fre[i] += 1;
			}
		}
	}	
	sum=0;
	for(i=1;i<Nbins;i++){
		sum += fre[i];
	}
	flag=true;
	flag1=true;
	if(normalize){
	for(i=1;i<Nbins;i++){
		if(i==1){
		pm[i] = (float) fre[i]/sum;
		Pm[i] =  pm[i];
		}
		else{
			pm[i] = (float) fre[i]/sum;
			Pm[i] = Pm[i-1] + pm[i];
		}
		if(Pm[i]>=0.05&&flag) {
			fprintf(fp7, "%f %f %f\n",
			(m3bins[i+1]+m3bins[i])/2.0,pm[i],Pm[i]); 
			flag=false;
		}
		if(Pm[i]>=0.95&&flag1){
			fprintf(fp7, "%f %f %f\n",(m3bins[i+1]+m3bins[i])/2.0,pm[i],Pm[i]);
			flag1=false;
		}
		fprintf(fp3,"%7.4f %6.4f %6.4f\n", (m3bins[i+1]+m3bins[i])/2.0,pm[i],Pm[i]);
	}
}
	else{
		for(i=1;i<Nbins;i++){
		if(i==1){
		pm[i] = (float) fre[i]/sum;
		Pm[i] =  pm[i];
		}
		else{
			pm[i] = (float) fre[i]/sum;
			Pm[i] = Pm[i-1] + pm[i];
		}
		if(Pm[i]>=0.05&&flag) {
			fprintf(fp7, "%f %f %f\n",
			(m3bins[i+1]+m3bins[i])/2.0,pm[i],Pm[i]); 
			flag=false;
		}
		if(Pm[i]>=0.95&&flag1){
			fprintf(fp7, "%f %f %f\n",(m3bins[i+1]+m3bins[i])/2.0,pm[i],Pm[i]);
			flag1=false;
		}
		fprintf(fp3,"%7.4f %6.4f %6.4f\n", (m3bins[i+1]+m3bins[i])/2.0,pm[i]*sum,Pm[i]*sum);
	}
}
	fclose(fp3);
	maxP3=max_vector(pm,Nbins-1);
	if(normalize) {
		fprintf(fp8,"%f %f\n",(maxP3.max_loc-1)*Dm3bins+(m3bins[2]+m3bins[1])/2.0,
		maxP3.max_value);
		fprintf(fp9,"%d\n",(int)(maxP3.max_value)/4);
	}
	else {
		fprintf(fp8,"%f %f\n",(maxP3.max_loc-1)*Dm3bins+(m3bins[2]+m3bins[1])/2.0,
		maxP3.max_value*sum);
		fprintf(fp9,"%d\n",(int)(maxP3.max_value*sum)/4);
	}
		for(i=1;i<Nbins;i++){
		fre[i] = 0;
		for(j=1;j<=count;j++){
			if(result[4][j]<m4bins[i+1]&&result[4][j]>=m4bins[i]){
				fre[i] += 1;
			}
		}
	}	
	sum=0;
	for(i=1;i<Nbins;i++){
		sum += fre[i];
	}
	flag=true;
	flag1=true;
	if(normalize){
	for(i=1;i<Nbins;i++){
		if(i==1){
		pm[i] = (float) fre[i]/sum;
		Pm[i] =  pm[i];
		}
		else{
			pm[i] = (float) fre[i]/sum;
			Pm[i] = Pm[i-1] + pm[i];
		}
		if(Pm[i]>=0.05&&flag) {
			fprintf(fp7, "%f %f %f\n",
			(m4bins[i+1]+m4bins[i])/2.0,pm[i],Pm[i]); 
			flag=false;
		}
		if(Pm[i]>=0.95&&flag1){
			fprintf(fp7, "%f %f %f\n",(m4bins[i+1]+m4bins[i])/2.0,pm[i],Pm[i]);
			flag1=false;
		}
		fprintf(fp4,"%7.4f %6.4f %6.4f\n", (m4bins[i+1]+m4bins[i])/2.0,pm[i],Pm[i]);
	}
}
	else{
		for(i=1;i<Nbins;i++){
		if(i==1){
		pm[i] = (float) fre[i]/sum;
		Pm[i] =  pm[i];
		}
		else{
			pm[i] = (float) fre[i]/sum;
			Pm[i] = Pm[i-1] + pm[i];
		}
		if(Pm[i]>=0.05&&flag) {
			fprintf(fp7, "%f %f %f\n",
			(m4bins[i+1]+m4bins[i])/2.0,pm[i],Pm[i]); 
			flag=false;
		}
		if(Pm[i]>=0.95&&flag1){
			fprintf(fp7, "%f %f %f\n",(m4bins[i+1]+m4bins[i])/2.0,pm[i],Pm[i]);
			flag1=false;
		}
		fprintf(fp4,"%7.4f %6.4f %6.4f\n", (m4bins[i+1]+m4bins[i])/2.0,pm[i]*sum,Pm[i]*sum);
	}
}
	fclose(fp4);
	maxP4=max_vector(pm,Nbins-1);
	if(normalize) {
		fprintf(fp8,"%f %f\n",(maxP4.max_loc-1)*Dm4bins+(m4bins[2]+m4bins[1])/2.0,
		maxP4.max_value);
		fprintf(fp9,"%d\n",(int)(maxP4.max_value)/4);
	}
	else {
		fprintf(fp8,"%f %f\n",(maxP4.max_loc-1)*Dm4bins+(m4bins[2]+m4bins[1])/2.0,
		maxP4.max_value*sum);
		fprintf(fp9,"%d\n",(int)(maxP2.max_value*sum)/4);
	}
	fclose(fp9);
	fclose(fp8);
	fclose(fp7);
	fprintf(fp5, "%f %f %f %f",Dm1bins,Dm2bins,Dm3bins,Dm4bins);
	fclose(fp5);
	fprintf(fp6, "%f %f %f %f %f %f %f %f",min1.min_value,max1.max_value
		,min2.min_value,max2.max_value,min3.min_value,max3.max_value,min4.min_value,max4.max_value);
	fclose(fp6);
	free_matrix(result,1,5,1,MAX_Array);
	free_vector(m1bins,1,Nbins);
	free_vector(m2bins,1,Nbins);
	free_vector(m3bins,1,Nbins);	
	free_vector(m4bins,1,Nbins);
	free_ivector(fre,1,Nbins);
	free_vector(pm,1,Nbins);
	free_vector(Pm,1,Nbins);
	return 0;
}
