#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

//An akima-spline struct is implemented to contain the data-points and the coefficients of the spline (see equation (1.30)).
typedef struct {gsl_vector *x, *y , *b, *c, *d;} akima_spline;

//The akima_spline_alloc function takes the two vectors with data-points and defines values for the b, c, and d-vectors in the struct.
akima_spline* akima_spline_alloc(gsl_vector *x, gsl_vector *y){
	akima_spline * s = (akima_spline*)malloc(sizeof(akima_spline)); // Memory is allocated to the akima_spline struct. 
	int n = x->size; // the number of data points is defined as n.
	
	//Memory is allocated to all vectors in the struct. The x and y vectors (data points) are defined as the x and y vectors in the struct.
	//The vectors b, c, and d contains the coefficients of equation (1.30) and they are defined in equation (1.32).
	s->x=gsl_vector_alloc(n);
	s->y=gsl_vector_alloc(n);
	s->b=gsl_vector_alloc(n);
	s->c=gsl_vector_alloc(n-1);
	s->d=gsl_vector_alloc(n-1);
	s->x=x; s->y=y;
	
	//Memory is allocated to the vectors h and p
	gsl_vector *h = gsl_vector_alloc(n-1);
	gsl_vector *p = gsl_vector_alloc(n-1);
	
	//The vector h is determined. h = Δx_i
	for(int i = 0; i < n-1; i++){
		gsl_vector_set(h,i,gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
		assert(gsl_vector_get(h,i)>0);
	}
	//Vector p is determined (defined just below equation (1.32)). p = Δy_i/Δx_i.
	for(int i = 0; i < n-1; i++){
		gsl_vector_set(p,i,(gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/gsl_vector_get(h,i));
	}
	
	//The first two and the last two points of vector b = s' is defined by equation (1.36):
	gsl_vector_set(s->b, 0, gsl_vector_get(p,0));
	gsl_vector_set(s->b, 1, (gsl_vector_get(p,0)+gsl_vector_get(p,1))/2);
	gsl_vector_set(s->b, n-1, gsl_vector_get(p,n-2));
	gsl_vector_set(s->b, n-2, (gsl_vector_get(p,n-2)+gsl_vector_get(p,n-3))/2);
	
	//The rest of the s' vector (b) is determined by equation (1.33) and (1.34), depending on the relation between w1 and w2. 
	//w1 and w2 are determined by equation (1.35).
	for(int i = 2; i < n-2; i++){
		double w1 = fabs(gsl_vector_get(p,i+1)-gsl_vector_get(p,i));
		double w2 = fabs(gsl_vector_get(p,i-1)-gsl_vector_get(p,i-2));
		if(w1+w2==0){gsl_vector_set(s->b, i, (gsl_vector_get(p,i-1)+gsl_vector_get(p,i))/2);}
		else gsl_vector_set(s->b, i, (w1*gsl_vector_get(p,i-1)+w2*gsl_vector_get(p,i))/(w1+w2));
	}
	
	//The c- and d-vectors are defined by equation (1.32).
	for(int i = 0; i < n-1; i++){
		gsl_vector_set(s->c, i, (3*gsl_vector_get(p,i)-2*gsl_vector_get(s->b, i)-gsl_vector_get(s->b, i+1))/gsl_vector_get(h,i));
		gsl_vector_set(s->d, i, (gsl_vector_get(s->b, i+1)+gsl_vector_get(s->b, i)-2*gsl_vector_get(p,i))/gsl_vector_get(h,i)/gsl_vector_get(h,i));
	}
	return s;
}

//The function akima_spline_eval calculates the spline-value at the point z.
double akima_spline_eval(akima_spline *s, double z){
	int n = s->x->size;
	assert(z>=gsl_vector_get(s->x, 0) && z <= gsl_vector_get(s->x, n-1)); //makes sure the z-value is inside the x-data-points values.
	int i = 0, j = n-1;
	while (j-i > 1){ //binary search for the interval for z
		int m = (i+j)/2;
		if (z>gsl_vector_get(s->x, m)) i = m;
		else j = m;
	}
	double h = z - gsl_vector_get(s->x, i); //distance from the z-point to the nearest x-point.
	//Calculates the interpolating spline:
	return gsl_vector_get(s->y, i) + h*(gsl_vector_get(s->b,i)+h*(gsl_vector_get(s->c,i)+h*gsl_vector_get(s->d,i)));
}

//The function akima_spline_free free the allocated memory.
void akima_spline_free(akima_spline *s){
	gsl_vector_free(s->x); gsl_vector_free(s->y); gsl_vector_free(s->b); gsl_vector_free(s->c); gsl_vector_free(s->d); free(s);
}