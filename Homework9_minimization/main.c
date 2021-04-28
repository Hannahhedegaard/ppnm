#include <math.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++)printf("%10g \n",gsl_vector_get(v,i));
//	printf("\n");
}

void qnewton(
	double f(gsl_vector* x), /* objective function */
	gsl_vector* x, /* on input: starting point, on exit: approximation to root */
	double eps /* accuracy goal, on exit |gradient| should be <eps */
	);

double f1(gsl_vector* x){
	double xx = gsl_vector_get(x,0);
    return pow((xx+1),2)-3;
}

double fRose(gsl_vector* x){
	double xx = gsl_vector_get(x,0);
	double yy = gsl_vector_get(x,1);
	return (pow(1-xx,2)+100*pow((yy-pow(xx,2)),2));
}


int main(){
	int n = 1;
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector_set(x,0,0);
	double eps = 1e-6;
	qnewton(f1, x, eps);
	vector_print("x =", x);
	
	//n = 2;
	//gsl_vector* x2 = gsl_vector_alloc(n);
	//gsl_vector_set(x2,0,1.1);
	//gsl_vector_set(x2,1,1.1);
	//qnewton(fRose, x2, eps);
	//vector_print("x_rose =", x2);
	
	
	return 0;
}