#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_blas.h>

void f1(gsl_vector* x,gsl_vector* fx){
    gsl_vector_set(fx,0,pow(gsl_vector_get(x,0)+1,2));        
}

void f2(gsl_vector* x,gsl_vector* fx){
    double vx = gsl_vector_get(x,0); 
    double vy = gsl_vector_get(x,1); 
    gsl_vector_set(fx,0,400*vx*vx*vx+(2-400*vy)*vx-2);        
    gsl_vector_set(fx,1,200*vy-200*vx*vx);      
}

double e; 
void f3(double r, gsl_vector* y, gsl_vector* dydx){
    gsl_vector_set(dydx,0,gsl_vector_get(y,1));  //
    gsl_vector_set(dydx,1,-2*gsl_vector_get(y,0)*(e+1/r));  //f'' = -2*f*(e+1/r)  
}

void driver(void f(double x,gsl_vector* y ,gsl_vector* dydx), double a, gsl_vector* ya,               
	double b, gsl_vector* yb, double h, double acc, double eps, char filename[]);

double rmax;
void g(gsl_vector* x, gsl_vector* fx){
    //rmax = 10 jf. opgaven
    e = gsl_vector_get(x, 0);
    double a = 1e-16, b = rmax, h = 1e-6, acc = 1e-5, ODEeps = 1e-4;
    gsl_vector* ya = gsl_vector_alloc(2);
    gsl_vector* yb = gsl_vector_alloc(2);
    //Initials values
    gsl_vector_set(ya,0,0); // Værdi af bølgefunktionen ved r=0, r-r^2 = 0 
    gsl_vector_set(ya,1,1); // Afledt ved r =0, f' = 1-2r => 1-0 = 1 
    driver(f3, a, ya, b, yb, h, acc, ODEeps, "hydrogenatom.txt");
    gsl_vector_set(fx,0,gsl_vector_get(yb,0)); 
}

void g2(gsl_vector* x, gsl_vector* fx){
    //rmax = 10 jf. opgaven
    e = gsl_vector_get(x, 0);
    double a = 1e-16, b = rmax, h = 1e-6, acc = 1e-5, ODEeps = 1e-4;
    gsl_vector* ya = gsl_vector_alloc(2);
    gsl_vector* yb = gsl_vector_alloc(2);
    //Initials values
    gsl_vector_set(ya,0,0); // Værdi af bølgefunktionen ved r=0, r-r^2 = 0 
    gsl_vector_set(ya,1,1); // Afledt ved r =0, f' = 1-2r => 1-0 = 1 
    driver(f3, a, ya, b, yb, h, acc, ODEeps, "hydrogenatom.txt");
    double k = pow(-2*e,0.5);
    gsl_vector_set(fx,0,gsl_vector_get(yb,0)-b*exp(-k*b)); 
}

void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps);
void vector_print(char s[], gsl_vector* v);

int main(){
    printf("Task A: Newton's method with numerical Jacobian and back-tracking linesearch\n");
    printf("Newton's method is implemented as the function newton in functions.c.\n");
    printf("The routine is checked by finding roots of the function: f(x)=(x+1)^2.\n");
    gsl_vector*x = gsl_vector_alloc(1);
    gsl_vector_set(x,0,5); 
    double eps = 0.0001;
    newton(f1,x,eps);
    vector_print("Root =", x);
    printf("The exact result is -1.\n");

    printf("Furthermore, the extremums of the Rosenbrock's valley function f(x,y)=(1-x)^2+100*(y-x^2)^2 is found with the implementation.\n");
    gsl_vector*x2 = gsl_vector_alloc(2);
    gsl_vector_set(x2,0,10);
    gsl_vector_set(x2,1,10);
    newton(f2,x2,eps);
    vector_print("Extremums =", x2);   
    printf("The exact result is (1,1).\n");

    gsl_vector_free(x);
    gsl_vector_free(x2);
    printf("---------------------------------------\n");
    printf("Task B: Bound states of hydrogen atom with shooting method for boundary value problems.\n");
    printf("The lowest root ε0 of the equation M(ε)=0 for rmax = 8 will now be calculated with the routines.\n");
    rmax = 8; 
    gsl_vector* x3 = gsl_vector_alloc(1);
    gsl_vector_set(x3,0,-1);
    eps = 1e-6; 
    newton(g,x3,eps);
    vector_print("ε0 =", x3);
    printf("The exact result is -1/2.\n");
    printf("The data from the calculation is written into hydrogenatom.txt.\n");
    printf("The plot of the resulting function is seen in Hydrogen.pyxplot.png.\n");
	
	gsl_vector_free(x3);
	
    printf("---------------------------------------\n");
    printf("Task C: Better boundary condition for the hydrogen atom problem.\n");
    printf("First, the boundary condition f(rmax) = 0 is used to investigate convergence of the solution as a function of rmax.\n");
    gsl_vector* x4 = gsl_vector_alloc(1);
    FILE * file = fopen("C_firstcondition.txt", "w"); 
    for(rmax = 2.5; rmax<15; rmax+=(double)1/2){
        gsl_vector_set(x4,0,-1);
        newton(g,x4,eps);
        fprintf(file,"%8.3g %.15g \n", rmax, gsl_vector_get(x4,0));
    }
    printf("The data from this first condition is written into C_firstcondition.txt.\n");
    printf("The data is plotted in C_firstcondition.pyxplot.png.\n");
    
    printf("\nNow, the boundary condition f(rmax->infinity)= r*exp(-k*r) is used instead, where k = sqrt(-2ε).\n");
    gsl_vector* x5 = gsl_vector_alloc(1);
    FILE * file2 = fopen("C_secondcondition.txt", "w"); 
    for(rmax = 2.5; rmax<15; rmax+=(double)1/2){
        gsl_vector_set(x5,0,-1);
        newton(g2,x5,eps);
        fprintf(file2,"%8.3g %.15g \n", rmax, gsl_vector_get(x5,0));
    }
    printf("The data from the second condition is written into C_secondcondition.txt.\n");
    printf("The data is plotted in C_secondcondition.pyxplot.png.\n");
    printf("For comparisson of the two see C1andC2.pyxplot.png.\n");
    printf("The plot shows that C_secondcondition converges much faster to the correct solution.\n");
	
    gsl_vector_free(x4);
    gsl_vector_free(x5);
	fclose(file);
    fclose(file2);
}