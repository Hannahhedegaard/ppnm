#include<assert.h>
#include<gsl/gsl_vector.h>
#include<math.h>

//Cubic spline taken from my Homework 1.
int binsearch(int n,gsl_vector *x, double z){/* locates the interval for z by bisection */ 
	assert(gsl_vector_get(x,0)<=z && z<=gsl_vector_get(x,n-1));
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>gsl_vector_get(x,mid)) i=mid; else j=mid;
		}
return i;
}

typedef struct {gsl_vector *x, *y, *b, *c, *d;} cubic_spline;

cubic_spline* cubic_spline_alloc(gsl_vector *x, gsl_vector* y){
	cubic_spline* s=(cubic_spline*)malloc(sizeof(cubic_spline));
	int n = x->size;

	s->b=gsl_vector_alloc(n);
	s->c=gsl_vector_alloc(n-1);
	s->d=gsl_vector_alloc(n-1);
	s->x=gsl_vector_alloc(n);
	s->y=gsl_vector_alloc(n);
	s->x=x; s->y=y;

	gsl_vector * h=gsl_vector_alloc(n-1);
	gsl_vector * p=gsl_vector_alloc(n-1);

	for(int i =0; i < n-1; i++)
		{
		gsl_vector_set(h,i,gsl_vector_get(x,i+1)-gsl_vector_get(x,i));
		assert(gsl_vector_get(h,i)>0);
		}

	for(int i =0; i < n-1; i++)
		{
		gsl_vector_set(p,i,(gsl_vector_get(y,i+1)-gsl_vector_get(y,i))/gsl_vector_get(h,i));
		}

	gsl_vector * D=gsl_vector_alloc(n);
	gsl_vector * Q=gsl_vector_alloc(n-1);
	gsl_vector * B=gsl_vector_alloc(n);

	gsl_vector_set(D,0,2);
	for(int i=0; i < n-2; i++)
		{
		gsl_vector_set(D,i+1,2*gsl_vector_get(h,i)/gsl_vector_get(h,i+1)+2);
		}

	gsl_vector_set(D,n-1,2);
	gsl_vector_set(Q,0,1);

	for(int i=0; i < n-2; i++)
		{
		gsl_vector_set(Q,i+1,gsl_vector_get(h,i)/gsl_vector_get(h,i+1));
		}

	for(int i = 0; i < n - 2; i++)
		{
		gsl_vector_set(B, i+1, 3*(gsl_vector_get(p,i)+gsl_vector_get(p,i+1)*gsl_vector_get(h,i)/gsl_vector_get(h,i+1)));
		}

	gsl_vector_set(B,0,3*gsl_vector_get(p,0));
	gsl_vector_set(B, n-1,3*gsl_vector_get(p,n-2));

	for (int i =1; i<n; i++)
		{
		gsl_vector_set(D,i,gsl_vector_get(D,i)-gsl_vector_get(Q,i-1)/gsl_vector_get(D,i-1));
		gsl_vector_set(B,i,gsl_vector_get(B,i)-gsl_vector_get(B,i-1)/gsl_vector_get(D,i-1));
		}

	gsl_vector_set(s->b, n-1, gsl_vector_get(B, n-1)/gsl_vector_get(D,n-1));

	for(int i=n-2; i>=0; i--)
		{
		gsl_vector_set(s->b, i, (gsl_vector_get(B,i)-gsl_vector_get(Q,i)*gsl_vector_get(s->b,i+1))/gsl_vector_get(D,i));
		}
	for(int i =0; i<n-1; i++)
		{
		gsl_vector_set(s->c, i, (-2*gsl_vector_get(s->b,i)-gsl_vector_get(s->b,i+1)+3*gsl_vector_get(p,i))/gsl_vector_get(h,i));
		gsl_vector_set(s->d, i, (gsl_vector_get(s->b,i)+gsl_vector_get(s->b,i+1)-2*gsl_vector_get(p,i))/gsl_vector_get(h,i)/gsl_vector_get(h,i));
		}

return s;
}

double cubic_spline_eval(cubic_spline*s, double z) {
	int index = binsearch(s->x->size, s->x, z);
	double h = z - gsl_vector_get(s->x,index);
	return gsl_vector_get(s->y,index)+h*(gsl_vector_get(s->b,index)+h*(gsl_vector_get(s->c,index)+h*(gsl_vector_get(s->d,index))));
}

void cubic_spline_free(cubic_spline *s){
	gsl_vector_free(s->b); gsl_vector_free(s->c); gsl_vector_free(s->d); free(s);
}