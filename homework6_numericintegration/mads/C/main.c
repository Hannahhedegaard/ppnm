#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

// struct for holding integration result and number of calls
typedef struct {
   double value;
   double error;
   int ncalls;
} integ;

// integration that also accepts infinite limits
void infintegrate(double f(double), double a, double b, double abs, double rel, integ* result);

// some functions to integrate
double f(double x) { // 0 to + inf
   return 2.0/(1 + x*x);
}
double g(double x) { // -inf to +inf
   return exp(-x*x);
}
double h(double x) { // -ing to 0
   return 3.0*pow(sin(x), 4)/pow(x, 4);
}


double ff(double x) { 
   return f(0 + x/(1-x))/pow(1 - x, 2);
}
double gg(double x) { 
   return g(x/(1 - x*x))*(1 + x*x)/pow(1 - x*x, 2);
}
double hh(double x) { 
   return h(x/(1 - x*x))*(1 + x*x)/pow(1 - x*x, 2);
}

// same functions but with dummy void* params to use with gsl integration routines
double fgsl(double x, void* params) { 
   return f(x);
}
double ggsl(double x, void* params) { 
   return g(x);
}
double hgsl(double x, void* params) { 
   return h(x);
}



// do stuff
int main() {


   // accurary goals
   double abs;
   double rel;

   // set-up for gsl integration
   int limit = 999; gsl_integration_workspace* w; w = gsl_integration_workspace_alloc(limit);
   double result, error;

   // integrate and print
   abs = 1e-6; rel = 1e-6; 
   integ Fstd; infintegrate(f, 0, INFINITY, abs, rel, &Fstd);
   //integ Fstd; infintegrate(ff, 0, 1, abs, rel, &Fstd);
   gsl_function Fgsl; Fgsl.function = &fgsl;
   gsl_integration_qagiu(&Fgsl, 0, abs, rel, limit, w, &result, &error);
   printf("--------------------------------------\n");
   printf("Integral of 2/(1 + x*x) from 0 to +inf\n");
   printf("--------------------------------------\n");
   printf("  Absolute accurary goal : %.1e\n" , abs);
   printf("  Relative accurary goal : %.1e\n" , rel);
   printf("  Exact value (pi)       : %.10f\n", M_PI);
   printf("  Variable trans. quad.  : \n");
   printf("    Calc. value          : %.10f\n", Fstd.value);
   printf("    Absolute difference  : %.1e\n" , fabs(Fstd.value - M_PI));
   printf("    Error estimate       : %.1e\n" , Fstd.error);
   printf("    Number of calls      : %d\n"   , Fstd.ncalls);
   printf("  GSL (qagiu)            : \n");
   printf("    Calc. value          : %.10f\n", result);
   printf("    Absolute difference  : %.1e\n" , fabs(result - M_PI));
   printf("    Error estimate       : %.1e\n\n" , error);

   abs = 1e-6; rel = 1e-6; 
   integ Gstd; infintegrate(g, -INFINITY, INFINITY, abs, rel, &Gstd);
   gsl_function Ggsl; Ggsl.function = &ggsl;
   gsl_integration_qagi(&Ggsl, abs, rel, limit, w, &result, &error);
   printf("---------------------------------------\n");
   printf("Integral of exp(-x*x) from -inf to +inf\n");
   printf("---------------------------------------\n");
   printf("  Absolute accurary goal : %.1e\n" , abs);
   printf("  Relative accurary goal : %.1e\n" , rel);
   printf("  Exact value (sqrt(pi)) : %.10f\n", sqrt(M_PI));
   printf("  Variable trans. quad.  : \n");
   printf("    Calc. value          : %.10f\n", Gstd.value);
   printf("    Absolute difference  : %.1e\n" , fabs(Gstd.value - sqrt(M_PI)));
   printf("    Error estimate       : %.1e\n" , Gstd.error);
   printf("    Number of calls      : %d\n"   , Gstd.ncalls);
   printf("  GSL (qagi)             : \n");
   printf("    Calc. value          : %.10f\n", result);
   printf("    Absolute difference  : %.1e\n" , fabs(result - sqrt(M_PI)));
   printf("    Error estimate       : %.1e\n\n" , error);


   abs = 1e-6; rel = 1e-6; 
   integ Hstd; infintegrate(h, -INFINITY, 0, abs, rel, &Hstd);
   gsl_function Hgsl; Hgsl.function = &hgsl;
   gsl_integration_qagil(&Hgsl, 0, abs, rel, limit, w, &result, &error);
   printf("-------------------------------------------\n");
   printf("Integral of 3*sin(x)**4/x**4 from -inf to 0\n");
   printf("-------------------------------------------\n");
   printf("  Absolute accurary goal : %.1e\n" , abs);
   printf("  Relative accurary goal : %.1e\n" , rel);
   printf("  Exact value (pi)       : %.10f\n", M_PI);
   printf("  Variable trans. quad.  : \n");
   printf("    Calc. value          : %.10f\n", Hstd.value);
   printf("    Absolute difference  : %.1e\n" , fabs(Hstd.value - M_PI));
   printf("    Error estimate       : %.1e\n" , Hstd.error);
   printf("    Number of calls      : %d\n"   , Hstd.ncalls);
   printf("  GSL (qagil)            : \n");
   printf("    Calc. value          : %.10f\n", result);
   printf("    Absolute difference  : %.1e\n" , fabs(result - M_PI));
   printf("    Error estimate       : %.1e\n\n" , error);

   gsl_integration_workspace_free(w);
  
   return 0;
}


