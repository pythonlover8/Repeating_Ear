#ifndef _takeoff_
#define _takeoff_

#include<stdio.h>
#include<math.h>
#include<string.h>
#include "util_err.h"
#define MAX_Array 16000
#define Max_file 1000
#define pi 3.1415926

typedef struct source {
	float vp,ev_radius;
}source;

extern void takeoff(float evdp, float a_radius, char *file, char *filename);
extern source model_search(float evdp, float a_radius, char *file);

#endif 
