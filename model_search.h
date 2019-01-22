#ifndef _MODEL_SEARCH_H_
#define _MODEL_SEARCH_H_

typedef struct{
	float vp;
	float ev_radius;
}source;
typedef struct max
{
	float max_value;
	float max_loc;
}max;
typedef struct min
{
	float min_value;
	float min_loc;
}min;
#define pi 3.1415926
source model_search(float evdp, float a_radius, char *file);
float* linespace(float fl, float fr, int n);
float* range(float fl, float fr, float interval);
int** unirand(int min_num, int max_num, int N, int M);
float* getcol(float **mat,int col);
float* get_vector_mat(float **mat, int **index, int count, int col);
float* get_vector_vec(float *vec, int **index, int count);
max max_vector(float *v, int count);
min min_vector(float *v, int count);
void swap_float(float xp,float yp);
void swap_int(int xp,int yp);
void sort_float(float *v,int N);
void sort_int(int *v,int N);
float median(float *x,int N);
#endif /* _MODEL_SEARCH_H */