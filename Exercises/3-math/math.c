#include<math.h>
#include<complex.h>
#include<stdio.h>

int main(){
	printf("Exercise 1\n");
	int n = 5;
	double x = 0.5;
	complex z = csqrt(-2);
	complex im=csqrt(-1);
	complex zz = cexp(im*M_PI);
	complex z2 = cexp(im);
	complex z3 = cpow(im, M_E);
	complex z4 = cpow(im, im);
	printf("Gamma function:\n");
	printf("Gamma(5) = %g\n", gamma(n));
	printf("\n");
	printf("Bessel function:\n");
	printf("J1(0.5) = %g\n", j1(x));
	printf("\n");
	printf("Checking functions:\n");
	printf("sqrt(-2) = %g + I*%g\n", creal(z), cimag(z));
	printf("e^(i*pi) = %g + I*%g\n", creal(zz), cimag(zz));
	printf("e^(i) = %g + I*%g\n", creal(z2), cimag(z2));
	printf("i^(e) = %g + I*%g\n", creal(z3), cimag(z3));
	printf("i^(i) = %g + I*%g\n", creal(z4), cimag(z4));
	printf("\n");
	printf("Exercise 2\n");
	printf("Finding the number of significant digits of\ntype float, double, and long double,\n");
	printf("using 1/9 = 0.1111111111... as example\n");
	printf("\n");
	float x_float = 1.f/9;
	double x_double = 1./9;
	long double x_long_double = 1.L/9;
	printf("Precision of float: %.25g\n", x_float);
	printf("Precision of double: %.25lg\n", x_double);
	printf("Precision of long double: %.25Lg\n", x_long_double);

return 0;
}
