#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_integration.h>

double adapt(double f(double x), double a, double b, double acc, double eps);
double infintegrate(double f(double x), double a, double b, double acc, double eps);

int calls=0;
double fu(double x){calls++; return 2.0/(1.0 + x*x);};
double gsl_fu(double x, void * params){calls++; return 2.0/(1.0 + x*x);};
//double f2(double x, void * params){
//	double f2 = 4*sqrt(1-pow(x,2));
//	return f2;
//};

double fb(double x){calls++; return 4*sqrt(1-x*x);};

double g(double x){return fb(cos(x))*sin(x);}; 	
double fuu(double x){
			calls++;
			return 2.0/(1.0 + pow((0 + (1.0 - x)/x),2)/(x*x));
		}


double g_gsl(double x, void * params){
	double g_gsl = fb(cos(x))*sin(x);
	return g_gsl;
}; 	

int main(){
	// opgave A
	double a=0, b=1; 
	double acc=1e-6, eps=1e-6;
	double Q1 = infintegrate(fuu,a,b,acc,eps); 
	printf("Q1=%.25g calls=%d\n", Q1, calls);
	
	// Opgave B
	calls = 0;
	a = 0; b = 1;
	a = acos(b); // Laver variable om
	b = acos(a); // Bytter om på grænserne fordi dx = -sinθdθ (se billede + afl-beskrivelse) - gælder uanset hvad a og b er (ikke inf)
	printf("a = %g, b = %g\n", a, b);
	double Qcc = infintegrate(g,a,b,acc,eps); 
	printf("Qcc=%.25g calls=%d\n", Qcc, calls);
	size_t size = 10000;
	double result, abserr;
	//sammenlign med GSL
	calls = 0;
	gsl_function F;
    F.function=&g_gsl;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(size);
	gsl_integration_qags(&F, a, b, acc, eps, size, w, &result, &abserr);
	printf("gsl=%.25g calls=%d\n", result, calls);
	
	//Opgave C
	calls = 0;
	b = INFINITY;
	double Q_inf = infintegrate(fu,a,b,acc,eps);
	printf("Q_inf = %.25g calls=%d\n", Q_inf, calls);
	//sammenlign med GSL
	calls = 0;
    F.function=&gsl_fu;
	gsl_integration_qagiu(&F, a, acc, eps, size, w, &result, &abserr);
	printf("gsl=%.25g calls=%d\n", result, calls);
	
return 0;
}



