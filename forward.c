#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"nrutil.h"
#include"model_search.h"


/* Input: x, y, z, t0, evdp*/
int main(int argc,char **argv){
	int i,j;
	int sto,sto1,sa,sa1;
	float x,y,z,origin_time;
	float travel_time,evdp,a_radius=6371.1;
	float *takeoff,*takeoff1,*azi,*azi1;
	float **result,**result1,time_interval,time_interval1;
	float t0,az,max_time=-10000,min_time=10000;
	char filename[160];
	source sour;
	FILE *fp,*fp1,*fp2,*fp3,*fp4,*fp5;
	fp=fopen("data.dat","w");
	if(fp==NULL){
		nrerror("unable to open the file");
	}
	fp1=fopen("data1.dat","w");
	if(fp1==NULL){
		nrerror("unable to open the file");
	}
	fp2=fopen("result33.dat","r");
	if(fp2==NULL){
		nrerror("unable to open the file");
	}
	fp3=fopen("result33pP.dat","r");
	if(fp3==NULL){
		nrerror("unable to open the file");
	}
	fp4=fopen("maxmin.dat","w");
	if(fp4==NULL){
		nrerror("unable to open the file");
	}
	fp5=fopen("time_interval.dat","w");
	if(fp5==NULL){
		nrerror("unable to open the file");
	}
	result=matrix(1,MAX_Size,1,5);
	result1=matrix(1,MAX_Size,1,5);
	for(i=1;!feof(fp2);i++){
		fscanf(fp2,"%f %f %f %f %f", &result[i][1], &result[i][2], &result[i][3],
			 &result[i][4], &result[i][5]);
		if(result[i][3]==0) continue;
		if(result[i][3]<min_time){
			min_time=result[i][3];
		}
		if(result[i][3]>max_time){
			max_time=result[i][3];
		}
	}
	fclose(fp2);
	for(i=1;!feof(fp3);i++){
		fscanf(fp3,"%f %f %f %f %f", &result1[i][1], &result1[i][2], &result1[i][3],
			 &result1[i][4], &result1[i][5]);
		if(result1[i][3]==0) continue;
		if(result1[i][3]<min_time){
			min_time=result1[i][3];
		}
		if(result1[i][3]>max_time){
			max_time=result1[i][3];
		}
	}
	fclose(fp3);

	if(argc<6){
		nrerror("Too few arguments!");
	}
	x = atof(argv[1]);
	y = atof(argv[2]);
	z = atof(argv[3]);
	origin_time = atof(argv[4]);
	evdp = atof(argv[5]);
    strcpy(filename,"/Users/a123/Desktop/Ear_relo/iasp91.dat");
	sour = model_search(evdp,a_radius,filename);
	sto = 91;
	sa = 361;
	takeoff = vector(1,sto);
	takeoff1 = vector(1,sto-1);
	azi = vector(1,sa);
	takeoff = linespace(0.0,90.0,sto);
	takeoff1 = linespace(90.0,180.0,sto);
	azi = linespace(0.0,360.0,sa);
/*	write data.dat for the background values */
	for(i=1;i<=sto;i++){
		t0 = takeoff[i];
		for(j=1;j<=sa;j++){
			az = azi[j];
			travel_time = -x*cos(az*pi/180)*sin(t0*pi/180)
				-y*sin(az*pi/180)*sin(t0*pi/180)+z*cos(t0*pi/180);
			travel_time = travel_time/sour.vp + origin_time;
			if(travel_time<min_time){
				min_time=travel_time;
			}
			if(travel_time>max_time){
				max_time=travel_time;
			}
			fprintf(fp, "%5.1f %5.1f %13.8f\n", az,t0,travel_time);	
		}
	}
	printf("%f %f\n",min_time,max_time);
	fclose(fp);
/* write data1.dat */
	for(i=1;i<=sto;i++){
		t0 = takeoff1[i];
		for(j=1;j<=sa;j++){
			az = azi[j];
			travel_time = -x*cos(az*pi/180)*sin(t0*pi/180)
				-y*sin(az*pi/180)*sin(t0*pi/180)+z*cos(t0*pi/180);
			travel_time = travel_time/sour.vp + origin_time;
			if(travel_time<min_time){
				min_time=travel_time;
			}
			if(travel_time>max_time){
				max_time=travel_time;
			}
			fprintf(fp1, "%5.1f %5.1f %13.8f\n", az,180 - t0,travel_time);
		}
	}
	fclose(fp1);
	printf("%f %f\n",min_time,max_time);
	time_interval=(max_time-min_time)/1000;
	time_interval1=(max_time-min_time)/6;
	fprintf(fp4,"%f %f",min_time,max_time);
	fprintf(fp5, "%f %f",time_interval,time_interval1);
	fclose(fp4);
	fclose(fp5);

	free_vector(takeoff,1,sto);
	free_vector(takeoff1,1,sto);
	free_vector(azi,1,sa);
	free_matrix(result,1,MAX_Size,1,5);
	free_matrix(result1,1,MAX_Size,1,5);
	return 0;
}
	