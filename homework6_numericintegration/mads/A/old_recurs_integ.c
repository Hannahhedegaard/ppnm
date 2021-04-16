#include<stdio.h>
#include<math.h>


// second order adaptive, recursive integrator
double integrate(double f(double), double a, double b, double fa, double fb, double abs, double rel, int count) {
   count += 1;
   double h = (b - a)/2;
   double fmid = f((a+b)/2);
   double Q = (h/3)*(fa + 4*fmid + fb);
   double q = (h/4)*(fa + 2*fmid + fb);
   double err = fabs(Q-q);
   if (err < abs + rel*fabs(Q)) return Q;
   else return integrate(f, a, (a+b)/2, fa, fmid, abs/sqrt(2), rel, count) + integrate(f, (a+b)/2, b, fmid, fb, abs/sqrt(2), rel, count);
}
