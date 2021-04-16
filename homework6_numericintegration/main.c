#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_integration.h>

double adapt(double f(double), double a, double b, double acc, double eps);
double infintegrate(double f(double), double a, double b, double acc, double eps);

int calls=0;
double fu(double x){calls++; return 2.0/(1.0 + x*x);};
//double f2(double x, void * params){
//	double f2 = 4*sqrt(1-pow(x,2));
//	return f2;
//};

//double g(double x){calls++; return f(cos(x))*sin(x);}; 	

double g2(double x, void * params){
	double g2 = fu(cos(x))*sin(x);
	return g2;
}; 	

int main(){
	double a=0.0, b=INFINITY; 
	double acc=0.01, eps=0.01;
	double Q = infintegrate(fu,a,b,acc,eps);
	printf("Q=%.25g calls=%d\n", Q, calls);
	//calls = 0;
	//a = acos(b); // Laver variable om
	//b = acos(a);
	//double Qcc = adapt(g,a,b,acc,eps);
	//printf("Qcc=%.25g calls=%d\n", Qcc, calls);
	//size_t size = 10000;
	//double result, abserr;
	//calls = 0;
	//gsl_function F;
    //F.function=&g2;
	//gsl_integration_workspace * w = gsl_integration_workspace_alloc(size);
	//gsl_integration_qags(&F, a, b, acc, eps, size, w, &result, &abserr);
	//printf("gsl=%.25g calls=%d\n", result, calls);
return 0;
}



