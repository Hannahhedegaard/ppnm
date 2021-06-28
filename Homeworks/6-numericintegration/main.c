#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_integration.h>

double adapt(double f(double x), double a, double b, double acc, double eps);
double infintegrate(double f(double x), double a, double b, double acc, double eps);

int calls=0;
double f1(double x){calls++; return sqrt(x);};
double f2(double x){calls++; return 4*sqrt(1-x*x);};
double f3(double x){calls++; return 1/sqrt(x);};
double g_1(double x){return f2(cos(x))*sin(x);}; 
double g_2(double x){return f3(cos(x))*sin(x);}; 
double f4(double x){calls++; return 2.0/(1.0 + x*x);};
		
double gsl_f4(double x, void * params){calls++; return 2.0/(1.0 + x*x);};

double g_gsl(double x, void * params){
	double g_gsl = f2(cos(x))*sin(x);
	return g_gsl;
}; 	

int main(){
	printf("Task A: Recursive adaptive integrator.\n");
	printf("A recursive adaptive integrator is implemented in functions.c.\n");
	printf("The integrator estimates the integral of a given function f(x)\non a given interval [a,b] with the required absolute, acc, or relative, eps, accuracy goals.\n");
	double a=0, b=1; 
	double acc=0.001, eps=0.001;
	printf("The implementation is tested on the two following functions:\n\n");
	
	printf("Integration of sqrt(x) from 0 to 1.\n");
	printf("acc = %g, eps = %g.\n", acc, eps);
	double Q1 = infintegrate(f1,a,b,acc,eps); 
	printf("Calculated integral = %.10g\n", Q1); 
	printf("Exact integral = %.10g\n", 2.0/3); 
	printf("Number of calls = %d\n\n", calls);
	calls = 0;
	printf("Integration of 4*sqrt(1-x^2) from 0 to 1.\n");
	printf("acc = %g, eps = %g.\n", acc, eps);
	double Q2 = infintegrate(f2,a,b,acc,eps); 
	printf("Calculated integral = %.10g\n", Q2); 
	printf("Exact integral = %.10g\n", M_PI); 
	printf("Number of calls = %d\n\n", calls);
	calls = 0;
	printf("Integration of 1/sqrt(x) from 0 to 1.\n");
	printf("acc = %g, eps = %g.\n", acc, eps);
	double Q3 = infintegrate(f3,a,b,acc,eps); 
	printf("Calculated integral = %.10g\n", Q3); 
	printf("Exact integral = %d\n", 2); 
	printf("Number of calls = %d\n", calls);
	
	printf("---------------------------------------\n");
	printf("Task B: Open quadrature with Clenshaw–Curtis variable transformation.\n");
	calls = 0;
	a = 0; b = 1;
	a = acos(b); // Laver variable om
	b = acos(a); // Bytter om på grænserne fordi dx = -sinθdθ (se billede + afl-beskrivelse) - gælder uanset hvad a og b er (ikke inf)
	printf("The new implementation is used on the two last functions from Task A, to compare accuracy and number of calls.\n");
	printf("The two new limits are calculated as:\n");
	printf("a = %g, b = %g\n\n", a, b);
	
	printf("Integration of 4*sqrt(1-x^2) from 0 to 1 with Clenshaw–Curtis variable transformation.\n");
	printf("acc = %g, eps = %g.\n", acc, eps);
	double Q4 = infintegrate(g_1,a,b,acc,eps); 
	printf("Calculated integral = %.10g\n", Q4); 
	printf("Exact integral = %.10g\n", M_PI); 
	printf("Number of calls = %d\n\n", calls);
	calls = 0;
	printf("Integration of 1/sqrt(x) from 0 to 1 with Clenshaw–Curtis variable transformation.\n");
	printf("acc = %g, eps = %g.\n", acc, eps);
	double Q5 = infintegrate(g_2,a,b,acc,eps); 
	printf("Calculated integral = %.10g\n", Q5); 
	printf("Exact integral = %d\n", 2); 
	printf("Number of calls = %d\n", calls);
	printf("Comment: in both cases, the number of calls have been greatly improved.\nThe accuracy is improved for the first integration but not significantly for the last.\n");
	
	printf("\nThe integration of 4*sqrt(1-x^2) is compared with the GSL routine:\n");
	size_t size = 10000;
	double result, abserr;
	calls = 0;
	gsl_function F;
    F.function=&g_gsl;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(size);
	gsl_integration_qags(&F, a, b, acc, eps, size, w, &result, &abserr);
	printf("Calculated GSL = %.10g\n", result);
	printf("Exact integral = %.10g\n", M_PI); 
	printf("Number of calls = %d\n", calls);
	printf("---------------------------------------\n");
	printf("Task C: Infinite limits\n");
	printf("The integrator is generalized to accept infinite limits.\n");
	printf("The implementation is tested on the function:\n");
	calls = 0;
	b = INFINITY;
	printf("Integration of 2/(1+x^2) from 0 to infinity.\n");
	printf("acc = %g, eps = %g.\n", acc, eps);
	double Q6 = infintegrate(f4,a,b,acc,eps); 
	printf("Calculated integral = %.10g\n", Q6); 
	printf("Exact integral = %.10g\n", M_PI); 
	printf("Number of calls = %d\n\n", calls);
	
	printf("Comparison with GSL routine:\n");
	calls = 0;
    F.function=&gsl_f4;
	gsl_integration_qagiu(&F, a, acc, eps, size, w, &result, &abserr);
	printf("Calculated GSL = %.10g\n", result);
	printf("Exact integral = %.10g\n", M_PI); 
	printf("Number of calls = %d", calls);
	
return 0;
}



