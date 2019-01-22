#ifndef _GRID_SEARCH_H_
#define _GRID_SEARCH_H_
#include "model_search.h"
#include<stdbool.h>
typedef struct {
	float min_err;
	float min_x;
	float min_y;
	float min_z;
	float min_t;
} err;

err grid_search(float xl, float xr, float grid_size_x,
	 float yl, float yr, float grid_size_y, float zl, float zr, 
	 float grid_size_z, float tl, float tr, float grid_size_t, 
	 float **result, source sour,int count, bool OUTPUT);

float min(float *result,int count);
float max(float *result,int count);
#endif /*_GRID_SEARCH_H*/	