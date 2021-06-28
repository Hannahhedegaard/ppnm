#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>

void rkstep12(void f(double x,gsl_vector* y,gsl_vector* dydx), double x, gsl_vector* yx,    
	double h, gsl_vector* yh, gsl_vector* dy); 
    
void driver(void f(double x,gsl_vector* y ,gsl_vector* dydx), double a, gsl_vector* ya,               
	double b, gsl_vector* yb, double h, double acc, double eps, char filename[]);

void f1(double x, gsl_vector* y, gsl_vector* dydx){
    gsl_vector_set(dydx,0,gsl_vector_get(y,1)); 
    gsl_vector_set(dydx,1,-gsl_vector_get(y,0));
}

void fSIR(double x, gsl_vector* y, gsl_vector* dydx){
    int N = 1000; 
    int Tc = 4; 
    int Tr = 10; 
    gsl_vector_set(dydx,0,-gsl_vector_get(y,1)*gsl_vector_get(y,0)/(N*Tc));
    gsl_vector_set(dydx,1,gsl_vector_get(y,1)*gsl_vector_get(y,0)/(N*Tc)-gsl_vector_get(y,1)/Tr);
    gsl_vector_set(dydx,2,gsl_vector_get(y,1)/Tr);  
}

void fNewton(double x, gsl_vector* y, gsl_vector* dydx){
    gsl_vector_set(dydx,0,gsl_vector_get(y,6));
    gsl_vector_set(dydx,1,gsl_vector_get(y,7));
    gsl_vector_set(dydx,2,gsl_vector_get(y,8));
    gsl_vector_set(dydx,3,gsl_vector_get(y,9));
    gsl_vector_set(dydx,4,gsl_vector_get(y,10));
    gsl_vector_set(dydx,5,gsl_vector_get(y,11));
    
    double y0 = gsl_vector_get(y,0);
    double y1 = gsl_vector_get(y,1);
    double y2 = gsl_vector_get(y,2);
    double y3 = gsl_vector_get(y,3);
    double y4 = gsl_vector_get(y,4);
    double y5 = gsl_vector_get(y,5);
    
    double R1 = pow(pow(y2-y0,2)+pow(y3-y1,2),0.5);
    double R2 = pow(pow(y4-y0,2)+pow(y5-y1,2),0.5);
    double R3 = pow(pow(y4-y2,2)+pow(y5-y3,2),0.5);
      
    gsl_vector_set(dydx,6,-(y0-y2)/pow(R1,3)-(y0-y4)/pow(R2,3)); 
    gsl_vector_set(dydx,7,-(y1-y3)/pow(R1,3)-(y1-y5)/pow(R2,3)); 
    gsl_vector_set(dydx,8, -(y2-y4)/pow(R3,3)-(y2-y0)/pow(R1,3)); 
    gsl_vector_set(dydx,9, -(y3-y5)/pow(R3,3)-(y3-y1)/pow(R1,3)); 
    gsl_vector_set(dydx,10,-(y4-y0)/pow(R2,3)-(y4-y2)/pow(R3,3)); 
    gsl_vector_set(dydx,11,-(y5-y1)/pow(R2,3)-(y5-y3)/pow(R3,3)); 
}


int main() {
    printf("Task A: Embedded Runge-Kutta ODE integrator\n");
    printf("The Runge-Kutta stepper with XY=12 is implemented in functions.c.\n");
    printf("Furthermore, an adaptive-step-size driver routine is implemented in functions.c.\n");
    printf("The routines are checked with the example: u''=-u.\n"); 
     
    int n = 2; 
    gsl_vector* ya = gsl_vector_alloc(n); 
    gsl_vector* yb = gsl_vector_alloc(n); 
    gsl_vector_set(ya,0,0); //sin(0) = 0
    gsl_vector_set(ya,1,1); //sin'(0) = cos(0) = 1
    double a = 0, b = 7, h = 0.1, acc = 1e-2, eps = 1e-2;
    driver(f1, a, ya, b, yb, h, acc, eps, "uprime.txt"); 
    
    //SIR
    printf("As another example, the SIR model of epidemic development is considered.\n");
    n = 3; 
    gsl_vector* yaSIR = gsl_vector_alloc(n); 
    gsl_vector* ybSIR = gsl_vector_alloc(n); 
    gsl_vector_set(yaSIR,0,900); 
    gsl_vector_set(yaSIR,1,10); 
    gsl_vector_set(yaSIR,2,0); 
    a = 0, b = 150, h = 0.1, acc = 1e-3, eps = 1e-3;
    driver(fSIR, a, yaSIR, b, ybSIR, h, acc, eps, "SIR.txt"); 
	
	printf("\nTask B: Store the path.\n");
    printf("The path from the uprime example is written in the file uprime.txt and plotted in uprime.pyxplot.png.\n");
	printf("The path from the SIR example is written into SIR.txt and plotted in SIR.pyxplot.png.\n"); 
	
    //Newton
    printf("\nTask C: Newtonian gravitational three-body problem.\n");
    printf("The routines are now used on Newtonian gravitational three-body problem.\n");
    printf("The path is written in Newton.txt and plotted in Newton.pyxplot.png.\n"); 
    printf("Furthermore, figure 8 is reproduced from this solution and is plotted in Figure8_Newton.pyxplot.png.\n"); 

    n = 12; 
    gsl_vector* yaNewton = gsl_vector_alloc(n); 
    gsl_vector* ybNewton = gsl_vector_alloc(n); 
    gsl_vector_set(yaNewton,0, 0.97000436); 
    gsl_vector_set(yaNewton,1, -0.24308753); 
    gsl_vector_set(yaNewton,2, -0.97000436); 
    gsl_vector_set(yaNewton,3, 0.24308753); 
    gsl_vector_set(yaNewton,4, 0); 
    gsl_vector_set(yaNewton,5, 0); 
    gsl_vector_set(yaNewton,6, 0.93240737/2); 
    gsl_vector_set(yaNewton,7, 0.86473146/2); 
    gsl_vector_set(yaNewton,8, 0.93240737/2); 
    gsl_vector_set(yaNewton,9,0.86473146/2); 
    gsl_vector_set(yaNewton,10, -0.93240737); 
    gsl_vector_set(yaNewton,11, -0.86473146); 
    a = 0, b = 10, h = 0.01, acc = 1e-3, eps = 1e-3;
    driver(fNewton, a, yaNewton, b, ybNewton, h, acc, eps, "Newton.txt"); 
    
return 0; 
} 