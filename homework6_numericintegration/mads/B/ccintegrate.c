#include<stdio.h>
#include<math.h>
#include<assert.h>

// adaptive, recursive quadrature (open with four points)
double quad(double f(double), double a, double b, double f2, double f3, double abs, double rel, int nrecur);

// struct for holding integration result and number of calls
typedef struct {
   double value;
   int ncalls;
} integ;


// integration wrapper with Clenshaw-Curtis transformation (uses gcc nested functions!)
void ccintegrate(double f(double), double a, double b, double abs, double rel, integ* result) {
   
   // counter variables 
   int ncalls = 0;
   int nrecur = 0;

   // do variable tranformation such that [a, b] -> [-1, 1]
   double ftrans(double x) {
      return f((b-a)*x/2 + (b+a)/2)*(b-a)/2;
   }

   // do clenshaw-curtis transformation
   double fclenshaw(double x) {
      ncalls += 1;
      return ftrans(cos(x))*sin(x);
   }
   
   // calculate first points
   double f2 = fclenshaw(2*M_PI/6);
   double f3 = fclenshaw(4*M_PI/6);

   // calculate integral
   double value = quad(fclenshaw, 0, M_PI, f2, f3, abs, rel, nrecur);

   // write data to result struct
   result->value = value; 
   result->ncalls = ncalls;

   // return integral and number of calls
   //printf("Number of calls: %d\n", result->ncalls);
   //printf("Result: %g\n", result->value);

}



