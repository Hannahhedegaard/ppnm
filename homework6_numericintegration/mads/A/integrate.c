#include<stdio.h>
#include<math.h>
#include<assert.h>

// struct for holding integration result and number of calls
typedef struct {
   double value;
   int ncalls;
} integ;

// adaptive, recursive quadrature (open with four points)
double quad(double f(double), double a, double b, double f2, double f3, double abs, double rel, int nrecur) {
   assert(nrecur < 1e5);
   double f1 = f(a + 1*(b-a)/6);
   double f4 = f(a + 5*(b-a)/6);
   double Q = (b-a)*(2*f1 + f2 + f3 + 2*f4)/6;
   double q = (b-a)*(  f1 + f2 + f3 +   f4)/4;
   double err = fabs(Q-q);
   if (err < abs + rel*fabs(Q)) return Q;
   else return quad(f, a, (a+b)/2, f1, f2, abs/sqrt(2), rel, nrecur + 1) + quad(f, (a+b)/2, b, f3, f4, abs/sqrt(2), rel, nrecur + 1);
}

// integration wrapper
void integrate(double f(double), double a, double b, double abs, double rel, integ* result) {

   // counter variables
   int nrecur = 0;
   int ncalls = 0;

   // stupid function to easily count number of calls (uses gcc nested functions!)
   double fstupid(double x) {
      ncalls += 1;
      return f(x);
   }

   // calculate first points
   double f2 = fstupid(a + 2*(b-a)/6);
   double f3 = fstupid(a + 4*(b-a)/6);

   // calculate integral
   double value = quad(fstupid, a, b, f2, f3, abs, rel, nrecur);
   
   // write data to result struct
   result->value = value; 
   result->ncalls = ncalls;

}



