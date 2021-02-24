#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

double f (double t, void* params) {
	double f = 2/sqrt(M_PI)*exp(-pow(t,2));
	return f;
}

double myfun(double z){
	gsl_function F;
	F.function = &f;
	int limit = 999;
	gsl_integration_workspace* w;
	w = gsl_integration_workspace_alloc (limit);
	double a=0, epsabs=1e-6, epsrel=1e-6, result, error;
	gsl_integration_qags(&F, a, z, epsabs, epsrel, limit, w, &result, &error);
	gsl_integration_workspace_free(w);
	return result;
}

int main(){
	for(double z=-3.0;z<=3.0;z+=1.0/8)
		printf("%10g %10g\n", z, myfun(z));
	
return 0;
}
