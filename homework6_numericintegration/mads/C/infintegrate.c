#include<stdio.h>
#include<math.h>
#include<assert.h>

// adaptive, recursive quadrature (open with four points)
double qquad(double f(double), double a, double b, double f2, double f3, double abs, double rel, int nrecur, double* error) {
   assert(nrecur < 1e5);
   double f1 = f(a + 1*(b-a)/6);
   double f4 = f(a + 5*(b-a)/6);
   double Q = (b-a)*(2*f1 + f2 + f3 + 2*f4)/6;
   double q = (b-a)*(  f1 + f2 + f3 +   f4)/4;
   double err = fabs(Q-q);
   if (err < abs + rel*fabs(Q)) {
      *error = err;
      return Q;
   } else {
      return qquad(f, a, (a+b)/2, f1, f2, abs/sqrt(2), rel, nrecur + 1, error) + qquad(f, (a+b)/2, b, f3, f4, abs/sqrt(2), rel, nrecur + 1, error);
   }
}

// struct for holding integration result and number of calls
typedef struct {
   double value;
   double error;
   int ncalls;
} integ;


// integration wrapper that accepts infinite limits via variable substitution (uses gcc nested functions!)
void infintegrate(double f(double), double a, double b, double abs, double rel, integ* result) {
   
   // initialize variables
   int ncalls = 0;
   int nrecur = 0;
   double value, error, A, B, f2, f3;

   // do variable tranformation if necessary 
   if (isinf(a) && isinf(b)) {

      double ftrans(double x) { ncalls++; return (f((1 - x)/x) + f(-(1 - x)/x))/(x*x); }
      A = 0; B = 1;   
      f2 = ftrans(A + 2.0*(B-A)/6.0);
      f3 = ftrans(A + 4.0*(B-A)/6.0);
      value = qquad(ftrans, A, B, f2, f3, abs, rel, nrecur, &error);

   } else if (isinf(b)) {

      double ftrans(double x) { ncalls++; return f(a + (1 - x)/x)/(x*x); }
      A = 0; B = 1;   
      f2 = ftrans(A + 2.0*(B-A)/6.0);
      f3 = ftrans(A + 4.0*(B-A)/6.0);
      value = qquad(ftrans, A, B, f2, f3, abs, rel, nrecur, &error);

   } else if (isinf(a)) {

      double ftrans(double x) { ncalls++; return f(b - (1 - x)/x)/(x*x); }
      A = 0; B = 1;   
      f2 = ftrans(A + 2.0*(B-A)/6.0);
      f3 = ftrans(A + 4.0*(B-A)/6.0);
      value = qquad(ftrans, A, B, f2, f3, abs, rel, nrecur, &error);

   } else { 

      double ftrans(double x) { ncalls++; return f(x); }
      A = a; B = b;   
      f2 = ftrans(A + 2.0*(B-A)/6.0);
      f3 = ftrans(A + 4.0*(B-A)/6.0);
      value = qquad(ftrans, A, B, f2, f3, abs, rel, nrecur, &error);

   }

   // write data to result struct
   result->value = value; 
   result->error = error; 
   result->ncalls = ncalls;

}



