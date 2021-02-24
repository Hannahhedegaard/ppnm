#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

double f (double tau, void* params) {
	double x = *(double*)params;
	double f = 1/M_PI*cos(tau-x*sin(tau));
	return f;
}

double myfun(double x){
	gsl_function F;
	F.function = &f;
	int limit = 999;
	gsl_integration_workspace* w;
	w = gsl_integration_workspace_alloc (limit);
	double a=0, b=M_PI, epsabs=1e-6, epsrel=1e-6, result, error;
	int key=6;
	gsl_integration_qag(&F, a, b, epsabs, epsrel, limit,key, w, &result, &error);
	gsl_integration_workspace_free(w);
	return result;
}

int main(){
	for(double z=1.0;z<=32;z+=0.1/8)
		printf("%10g %10g\n", z, myfun(z));
	
return 0;
}
