#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

// struct for holding integration result and number of calls
typedef struct {
   double value;
   int ncalls;
} integ;

// direct quadrature (i.e. no variable transformation) 
void integrate(double f(double), double a, double b, double abs, double rel, integ* result);

// integration using clenshaw-curtis transformation
void ccintegrate(double f(double), double a, double b, double abs, double rel, integ* result);

// some functions to integrate
double f(double x) { 
   return 1/sqrt(x);
}
double g(double x) { 
   return log(x)/sqrt(x);
}
double h(double x) { 
   return 4*sqrt(1 - x*x);
}

// same functions but with dummy void* params to use with gsl integration routines
double fgsl(double x, void* params) { 
   return 1/sqrt(x);
}
double ggsl(double x, void* params) { 
   return log(x)/sqrt(x);
}
double hgsl(double x, void* params) { 
   return 4*sqrt(1 - x*x);
}



// do stuff
int main() {

   // accurary goals
   double abs;
   double rel;

   // set-up for gsl integration
   int limit = 999; gsl_integration_workspace* w; w = gsl_integration_workspace_alloc(limit);
   int key = 2; double result, error;

   // integrate and print
   abs = 1e-4; rel = 1e-4; 
   integ Fstd;   integrate(f, 0, 1, abs, rel, &Fstd);
   integ Fcci; ccintegrate(f, 0, 1, abs, rel, &Fcci);
   gsl_function Fgsl; Fgsl.function = &fgsl;
   gsl_integration_qags(&Fgsl, 0, 1, abs, rel, limit, w, &result, &error);
   printf("---------------------------------\n");
   printf("Integral of 1/sqrt(x) from 0 to 1\n");
   printf("---------------------------------\n");
   printf("  Absolute accurary goal : %.1e\n" , abs);
   printf("  Relative accurary goal : %.1e\n" , rel);
   printf("  Standard quadrature    : \n");
   printf("    Value (should be 2)  : %.10f\n", Fstd.value);
   printf("    Number of calls      : %d\n"   , Fstd.ncalls);
   printf("  Clenshaw-Curtis        : \n");
   printf("    Value (should be 2)  : %.10f\n", Fcci.value);
   printf("    Number of calls      : %d\n"   , Fcci.ncalls);  
   printf("  GSL (qags)             : \n");
   printf("    Value (should be 2)  : %.10f\n\n", result);

   abs = 1e-4; rel = 1e-4; 
   integ Gstd;   integrate(g, 0, 1, abs, rel, &Gstd);
   integ Gcci; ccintegrate(g, 0, 1, abs, rel, &Gcci);
   gsl_function Ggsl; Ggsl.function = &ggsl;
   gsl_integration_qags(&Ggsl, 0, 1, abs, rel, limit, w, &result, &error);
   printf("-------------------------------------\n");
   printf("Integral of ln(x)/sqrt(x) from 0 to 1\n");
   printf("-------------------------------------\n");
   printf("  Absolute accurary goal : %.1e\n" , abs);
   printf("  Relative accurary goal : %.1e\n" , rel);
   printf("  Standard quadrature    : \n");
   printf("    Value (should be -4) : %.10f\n", Gstd.value);
   printf("    Number of calls      : %d\n"   , Gstd.ncalls);
   printf("  Clenshaw-Curtis        : \n");
   printf("    Value (should be -4) : %.10f\n", Gcci.value);
   printf("    Number of calls      : %d\n" , Gcci.ncalls); 
   printf("  GSL (qags)             : \n");
   printf("    Value (should be -4) : %.10f\n\n", result);


   // loop over tolerances to test how far we can get
   integ Hstd; integ Hcci;
   gsl_function Hgsl; Hgsl.function = &hgsl;
   double delta[] = {1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
   printf("------------------------------------------------------\n");
   printf("Integral of 4*sqrt(1 - x*x) from 0 to 1 (should be pi)\n");
   printf("------------------------------------------------------\n");
   printf("std  = open four point quadrature (no variable transformation)\n");
   printf("cc   = clenshaw-curtis transformation before using std\n");
   printf("qags = gsl qags routine\n");
   printf("D    = absolute difference wrt. exact value of pi\n");
   printf("N    = number of function calls\n");
   printf("%10s %10s %10s %10s %10s %10s %10s\n", "abs. acc.", "rel. acc.", "D(std)", "D(cc)", "D(qags)", "N(std)", "N(cc)");
   for (int i = 0; i < 6; ++i) {
      abs = delta[i];
      rel = delta[i];
      integrate(h, 0, 1, abs, rel, &Hstd);
      ccintegrate(h, 0, 1, abs, rel, &Hcci);
      gsl_integration_qag(&Hgsl, 0, 1, abs, rel, limit, key, w, &result, &error);
      printf("%10.1e %10.1e %10.1e %10.1e %10.1e %10d %10d\n", abs, rel, fabs(Hstd.value - M_PI), 
            fabs(Hcci.value - M_PI), fabs(result - M_PI), Hstd.ncalls, Hcci.ncalls);
   }

   gsl_integration_workspace_free(w);
   return 0;
}


