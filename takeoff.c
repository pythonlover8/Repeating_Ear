#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include<stdbool.h>
#include "model_search.h"
#include "nrutil.h"
#include "sac.h"

void write_sachdr(char *sacfile,float parameter){
	SACHEAD hd;
	float *trace;
	if((trace=read_sac(sacfile,&hd))==NULL){
		nrerror("unable to open the sacfile");
	}
	hd.user0 = parameter;
	write_sac(sacfile,hd,trace);
	free(trace);
}

/*input evdepth, earth_radius, the file of rayparameters, and the list of sacfiles*/
int main(int argc,char **argv){
	int i;
	float *takeoff_angle=NULL,*rayP=NULL,mode;
	float a_radius,evdp;
	source sour;
	char filename[160],file1[160],m[160],saclist[100],sacfile[160];
	strcpy(filename,"/Users/a123/Desktop/Ear_relo/iasp91.dat");
	strcpy(m,argv[3]);
	strcpy(saclist,argv[4]);
	evdp = atof(argv[1]);
	a_radius = atof(argv[2]);
	sour = model_search(evdp,a_radius,filename);
	strcpy(file1,"takeoff_angle.dat");
	takeoff_angle = vector(1,MAX_Size);
	rayP = vector(1,MAX_Size);
	FILE *f1, *f2, *f3;
	if((f1=fopen(file1,"w"))==NULL||(f2=fopen(m,"r"))==NULL||(f3=fopen(saclist,"r"))==NULL){
		nrerror("unable to open the file");
	}
	printf("Mode: -0 represents upward phase \n -1 represents downward phase \n");
	fscanf(stdin,"%f",&mode);
	for(i=1;!feof(f2);i++){
		fscanf(f2,"%f",rayP+i);
		if(*(rayP+i)==0){
			continue;
		}
		fscanf(f3,"%s",sacfile);
		if(mode){
		*(takeoff_angle+i) = asin(180/pi*(rayP[i]*(sour.vp/sour.ev_radius)));
		*(takeoff_angle+i) = *(takeoff_angle+i)*180/pi;
		write_sachdr(sacfile,takeoff_angle[i]);
		}
		else{
			*(takeoff_angle+i) = asin(180/pi*(rayP[i]*(sour.vp/sour.ev_radius)));
			*(takeoff_angle+i) = 180 - takeoff_angle[i]*180/pi;
			write_sachdr(sacfile,takeoff_angle[i]);
		}
		fprintf(f1,"%8.4f\n",*(takeoff_angle+i));
		
	}
	free_vector(takeoff_angle,1,MAX_Size);
	free_vector(rayP,1,MAX_Size);
	fclose(f1);
	fclose(f2);
	fclose(f3);
	return 0;
}
