#include<stdio.h> 
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<math.h> 
#include<float.h>
#include "nrutil.h"
#include "model_search.h"
#include "grid_search.h"
err grid_search(float xl, float xr, float grid_size_x,
	 float yl, float yr, float grid_size_y, float zl, float zr, 
	 float grid_size_z, float tl, float tr, float grid_size_t, 
	 float **result, source sour,int count, bool OUTPUT){
	int i,ii,iii,iiii; 
	int nx,ny,nz,nt,j;
	float *x,*y,*z,*t,travel_time;
	float error=0.0,time_lag,az,t0;
	err err1;
	err1.min_err = 4000;
	err1.min_x = 1000;
	err1.min_y = 1000;
	err1.min_z = 1000;
	err1.min_t = 1000;

	nx = ceil((xr-xl)/grid_size_x);
	ny = ceil((yr-yl)/grid_size_y);
	nz = ceil((zr-zl)/grid_size_z);
	nt = ceil((tr-tl)/grid_size_t);
	x = vector(0,nx-1);
	y = vector(0,ny-1);
	z = vector(0,nz-1);
	t = vector(0,nt-1);
	x = range(xl,xr,grid_size_x);
	y = range(yl,yr,grid_size_y);
	z = range(zl,zr,grid_size_z);
	t = range(tl,tr,grid_size_t);
	
	
	for(i=0;i<nx;i++){
		for(ii=0;ii<ny;ii++){
			for(iii=0;iii<nz;iii++){
				for(iiii=0;iiii<nt;iiii++){
					error = 0.0;
					for(j=0;j<count;j++){
						az = result[j][1];
						t0 = result[j][4];
						time_lag = result[j][2];
						travel_time = -x[i]*cos(az*pi/180)*sin(t0*pi/180)
							-y[ii]*sin(az*pi/180)*sin(t0*pi/180)+z[iii]*cos(t0*pi/180);
						travel_time = travel_time/sour.vp + t[iiii];
						error += fabs(travel_time-time_lag);
					}
					error /= count;
					if(error<err1.min_err){
						err1.min_err = error;
						err1.min_x = x[i];
						err1.min_y = y[ii];
						err1.min_z = z[iii];
						err1.min_t = t[iiii];
					}
				}
			}
		}
	}
	
	free_vector(x,0,MAX_Array-1);
	free_vector(y,0,MAX_Array-1);
	free_vector(z,0,MAX_Array-1);
	free_vector(t,0,MAX_Array-1);
	return err1;
}
float min(float *result,int count){
	int i;
	float min_value=100000;
	for(i=0;i<count;i++){
		if(result[i]<min_value){
			min_value = result[i];
		}
	}
	return min_value;
}
float max(float *result,int count){
	int i;
	float max_value=-100000;
	for(i=0;i<count;i++){
		if(result[i]>max_value){
			max_value = result[i];
		}
	}
	return max_value;
}
