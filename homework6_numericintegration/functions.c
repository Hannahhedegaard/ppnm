#include <math.h>
#include <stdio.h>
#include <assert.h>

double global_error;
int calls;

double integrate(double f(double x), double a, double b, double acc, double eps, double f2, double f3, int nrec)
{
	assert(nrec<1e6);
	printf("test = %g\n", f(1.0));
	//double f1=f(a+(b-a)/6.0), f4=f(a+5.0*(b-a)/6.0);
	//double Q=(2.0*f1+f2+f3+2.0*f4)/6.0*(b-a), q=(f1+f4+f2+f3)/4.0*(b-a);
	//double tolerance=acc+eps*fabs(Q), error=fabs(Q-q);
	//if (error < tolerance) {
	//	global_error += pow(error,2);
	//	return Q;
	//}
	//else {return integrate(f,a,(a+b)/2.0,acc/sqrt(2),eps, f1, f2, nrec+1)+
	//			integrate(f,(a+b)/2.0,b,acc/sqrt(2),eps, f3, f4, nrec+1);
	//}
	return 13;
}

double adapt(double f(double), double a, double b, double acc, double eps){
	global_error = 0.0;
	//double f2=f(a+2.0*(b-a)/6.0), f3=f(a+4.0*(b-a)/6.0); int nrec=0;
	//printf("f2 = %g\n", f2);
	//double estimate = integrate(f, a, b, acc, eps, f2, f3, nrec);
	printf("error = %g\n", sqrt(global_error));
	return 0;//estimate;
}

double infintegrate(double f(double), double a, double b, double acc, double eps){
	if (isinf(a)==-1 && isinf(b)==1){
		double g(double x){
			calls++;
			return (f((1-x)/x) + f(-(1.0-x)/x))*1.0/pow(x,2);
		}
		
		return adapt(g, 0.0, 1.0, acc, eps);
	}
	else if (isinf(b)==1){
		double g(double x){
			calls++;
			return f(a + (1.0 - x)/x)/(x*x);
		}
		double A = 0.0, B = 1.0;
		double f2=g(A+2.0*(B-A)/6.0), f3=g(A+4.0*(B-A)/6.0); int nrec=0;
		printf("f2=%g, f3=%g\n", f2,f3);
		printf("hej\n");
		double Result = integrate(g, A, B, acc, eps, f2, f3, nrec);
		printf("result = %g\n", Result);
		return 0;
	}
	else if (isinf(a)==-1){
		double g(double x){
			calls++;
			return f(b - (1.0-x)/x)*1.0/pow(x,2);
		}
		return adapt(g, 0.0, 1.0, acc, eps);
	}
return adapt(f, a, b, acc, eps);
}