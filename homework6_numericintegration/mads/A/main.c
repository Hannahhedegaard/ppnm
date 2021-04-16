#include<stdio.h>
#include<math.h>

// struct for holding integration result and number of calls
typedef struct {
   double value;
   int ncalls;
} integ;

// adaptive, recursive integrator
void integrate(double f(double), double a, double b, double abs, double rel, integ* result);

// some function we would like to integrate
double f(double x) { 
   return sqrt(x);
}

// another function
double g(double x) {
   return 4*sqrt(1 - x*x);
}

// yet another function
double h(double x) {
   return sin(x);
}

int main() {

   // accuracy goals
   double abs;
   double rel;

   // integrate and print result
   abs = 1e-4; rel = 1e-4; 
   integ F; integrate(f, 0, 1, abs, rel, &F);
   printf("-------------------------------\n");
   printf("Integral of sqrt(x) from 0 to 1\n");
   printf("-------------------------------\n");
   printf("  Absolute accurary goal : %.1e\n"   , abs);
   printf("  Relative accurary goal : %.1e\n"   , rel);
   printf("  Value (should be 2/3)  : %.10f\n"  , F.value);
   printf("  Number of calls        : %d\n\n"     , F.ncalls);

   abs = 1e-4; rel = 1e-4; 
   integ G; integrate(g, 0, 1, abs, rel, &G);
   printf("---------------------------------------\n");
   printf("Integral of 4*sqrt(1 - x*x) from 0 to 1\n");
   printf("---------------------------------------\n");
   printf("  Absolute accurary goal : %.1e\n"   , abs);
   printf("  Relative accurary goal : %.1e\n"   , rel);
   printf("  Value (should be pi)   : %.10f\n"  , G.value);
   printf("  Number of calls        : %d\n\n"     , G.ncalls);

   abs = 1e-4; rel = 1e-4; 
   integ H; integrate(h, 0, M_PI, abs, rel, &H);
   printf("------------------------------\n");
   printf("Integral of sin(x) from 0 to 1\n");
   printf("------------------------------\n");
   printf("  Absolute accurary goal : %.1e\n"   , abs);
   printf("  Relative accurary goal : %.1e\n"   , rel);
   printf("  Value (should be 2)    : %.10f\n"  , H.value);
   printf("  Number of calls        : %d\n\n"     , H.ncalls);


   //
   //
   //printf("Integral of sqrt(x) from 0 to 1 (should be 2/3)       : %g\n", integrate(f, 0, 1, abs, rel));
   //printf("Integral of 4*sqrt(1 - x*x) from 0 to 1 (should be pi): %g\n", integrate(g, 0, 1, abs, rel));
   //printf("Integral of sin(x) from 0 to pi (should be 2)         : %g\n", integrate(h, 0, M_PI, abs, rel));

   return 0; 
}
