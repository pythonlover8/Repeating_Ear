#include<stdio.h> 
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h> 
#include<float.h>
#include "nrutil.h"
#include "matrix.h"
#include "model_search.h"
#include "grid_search.h"
#define MAX_Array1 2e8
void grid_search(float xl, float xr, float grid_size_x,
	 float yl, float yr, float grid_size_y, float zl, float zr, 
	 float grid_size_z, float tl, float tr, float grid_size_t, 
	 float **result, source sour,int count){
	int i,ii,iii,iiii,M=4; 
	int nx,ny,nz,nt,j;
	float *x,*y,*z,*t,travel_time;
	float time_lag,az,t0;
	float error;
	FILE *fp=fopen("L1result.dat","w");
	err err1;
	/*err *err2=NULL;*/
	err1.err = 4000;
	err1.x = 1000;
	err1.y = 1000;
	err1.z = 1000;
	err1.t = 1000;

	nx = ceil((xr-xl)/grid_size_x);
	ny = ceil((yr-yl)/grid_size_y);
	nz = ceil((zr-zl)/grid_size_z);
	nt = ceil((tr-tl)/grid_size_t);
	x = vector(1,nx);
	y = vector(1,ny);
	z = vector(1,nz);
	t = vector(1,nt);
	x = range(xl,xr,grid_size_x);
	y = range(yl,yr,grid_size_y);
	z = range(zl,zr,grid_size_z);
	t = range(tl,tr,grid_size_t);
	for(i=1;i<=nx;i++){
		for(ii=1;ii<=ny;ii++){
			for(iii=1;iii<=nz;iii++){
				for(iiii=1;iiii<=nt;iiii++){
					error = 0.0;
					for(j=1;j<=count;j++){
						az = result[j][2];
						t0 = result[j][5];
						time_lag = result[j][3];
						travel_time = -x[i]*cos(az*pi/180)*sin(t0*pi/180)
							-y[ii]*sin(az*pi/180)*sin(t0*pi/180)+z[iii]*cos(t0*pi/180);
						travel_time = travel_time/sour.vp + t[iiii];
						error += fabs(travel_time-time_lag);

					}
					error = error/(count-M);
/*					err2[k].err=error;
					err2[k].x=x[i];
					err2[k].y=y[ii];
					err2[k].z=z[iii];
					err2[k].t=t[iiii];
*/
					if(error<err1.err){
						err1.err = error;
						err1.x = x[i];
						err1.y = y[ii];
						err1.z = z[iii];
						err1.t = t[iiii];
					}

				}
			}
		
		}
	}
	fprintf(fp, "%f %f %f %f %f\n",err1.x,err1.y,err1.z,
		err1.t,err1.err);
	fclose(fp);
	free_vector(x,1,MAX_Array);
	free_vector(y,1,MAX_Array);
	free_vector(z,1,MAX_Array);
	free_vector(t,1,MAX_Array);
/*	return err2; */
}
/*input: searching range and the corresponding grid-size ---\n 
     		min x, max x, and grid_size_x, min y, max y, and grid_size_y,\n
     		min z, max z, and grid_size_z, and min t, max t, and grid_size_t;\n
     		the file storing the results; sour_velocity.\n
     		Output: grdfile, and model_parameter.*/
int main(int argc, char **argv){
	int i,j,count;
	float xl,xr,yl,yr,zl,zr,tl,tr;
	float dx,dy,dz,dt;
	float **result,**result1;
	float evdp,a_radius=6371.1;
	char *file;
	FILE *fp;
	source sour;
	result=matrix(1,5,1,MAX_Array);
	result1=matrix(1,MAX_Array,1,5);
	printf("Input: searching range and the corresponding grid-size ---\n \
     	 	min x, max x, and grid_size_x, min y, max y, and grid_size_y,\n \
     	 	min z, max z, and grid_size_z, and min t, max t, and grid_size_t;\n");
	scanf("%f %f %f\n\
		%f %f %f\n\
		%f %f %f\n\
		%f %f %f",&xl,&xr,&dx,&yl,&yr,&dy,&zl,&zr,&dz,&tl,&tr,&dt);
	printf("Input: source depth\n");
	scanf("%f",&evdp);
	printf("running...\n");
	sour=model_search(evdp,a_radius,"/Users/a123/Desktop/Ear_relo/iasp91.dat");
	strcpy(file,"result3.dat");
	if((fp=fopen(file,"r"))==NULL){
		nrerror("Unable to open file");
	}
	for(i=1;!feof(fp);i++){
		fscanf(fp,"%f %f %f %f %f",&result[1][i],&result[2][i],
			&result[3][i],&result[4][i],&result[5][i]);

	}
	count=i-3;
	fclose(fp);
	transpose(result,result1,5,count);
	printf("running...\n");
	grid_search(xl,xr,dx,yl,yr,dy,zl,zr,dz,tl,tr,dt,result1,sour,count);
	free_matrix(result,1,5,1,MAX_Array);
	free_matrix(result,1,MAX_Array,1,5);
	return 0;
	}

