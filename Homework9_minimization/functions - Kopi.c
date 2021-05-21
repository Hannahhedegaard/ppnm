#include <math.h>
#include <float.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
//#define DBL_EPSILON 2.22045e-16

void GS_decomp(gsl_matrix* A, gsl_matrix* R);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void matrix_print(char s[], gsl_matrix* M);
void vector_print(char s[], gsl_vector* v);


void qnewton(
	double f(gsl_vector* x), /* objective function */
	gsl_vector* x, /* on input: starting point, on exit: approximation to root */
	double eps /* accuracy goal, on exit |gradient| should be <eps */
	){
	int n = x->size;
	int iter = 0;
	double lambda;
	double dfx;
	double fx;
	double fxs;
	double df;
	double fy;
	double uy;
	double sy;
	double gamma;
	double s_grad;
	double dx = sqrt(DBL_EPSILON);
	
	gsl_matrix* B = gsl_matrix_alloc(n,n);
	gsl_matrix_set_identity(B);
    gsl_vector* grad_F = gsl_vector_alloc(n);         
    gsl_vector* grad_Fs = gsl_vector_alloc(n);         
    gsl_vector* Dx = gsl_vector_alloc(n);    
    gsl_vector* y = gsl_vector_alloc(n);            
    gsl_vector* By = gsl_vector_alloc(n);        
    gsl_vector* u = gsl_vector_alloc(n);    
    
    while(1){
		iter += 1;
		if (iter>1e6) break;
        fx = f(x); 
        for (int j =0; j<n; j++) {
            gsl_vector_set(x,j,gsl_vector_get(x,j)+dx);
            dfx = f(x);
			df = dfx - fx;
            gsl_vector_set(grad_F,j,df/dx);             
			gsl_vector_set(x,j,gsl_vector_get(x,j)-dx);  
        }
		gsl_blas_dgemv(CblasNoTrans, -1, B, grad_F, 0, Dx);
    
        lambda = 2;
        while(1){
            lambda = lambda/2.0;
            gsl_vector_scale(Dx,lambda);
            gsl_vector_memcpy(y,x);
            gsl_vector_add(y,Dx);
            fy = f(y);
			gsl_blas_ddot(Dx, grad_F, &s_grad);
            if(fy<fx+1e-4*s_grad){
                break;
            }
			
			if (lambda < dx) {
				gsl_matrix_set_identity(B);
				break;
			}
			gsl_vector_scale(Dx,1.0/lambda);
        }
        
        if(gsl_blas_dnrm2(Dx) < dx*gsl_blas_dnrm2(x)) {
                break;
                }
		if(gsl_blas_dnrm2(grad_F) < eps) {
                break;
                }
	
	// Opdatering af B til B + δB
	gsl_vector_scale(Dx,lambda); // s-vektor bestemmes
	gsl_vector_add(Dx, x); 	// s <- s + x
	fxs = f(Dx); 			// funktionsværdien i (s+x)-vektoren
	for (int j =0; j<n; j++) {
		gsl_vector_set(Dx,j,gsl_vector_get(Dx,j)+dx);
        dfx = f(Dx); 		// funktionsværdi i (s+x) + dx
		df = dfx - fxs; 	// forskel i fkt.-værdi før og efter dx er lagt til
        gsl_vector_set(grad_Fs,j,df/dx); //bestemmer gradienten af s+x            
		gsl_vector_set(Dx,j,gsl_vector_get(Dx,j)-dx);   //går dx tilbage
    }
	
	// u og y-vektor bestemmes
	gsl_vector_sub(Dx, x); // Vektoren Dx er lig s. (s <- s - x)	
	gsl_vector_memcpy(u,Dx); // u = s
	gsl_vector_sub(grad_Fs, grad_F); // y-vektor bestemmes y = grad_Fs - grad_F
	gsl_blas_dgemv(CblasNoTrans, 1.0, B, grad_Fs, 0.0, By); // B*y -> By
	gsl_vector_sub(u, By); // u = s - By
	// prikprodukter uy og sy bestemmes
	
	gsl_blas_ddot(Dx, grad_Fs, &sy); // prikprodukt s^T*y = sy
	
	if (fabs(sy)>1e-12){ // laver B om hvis sy er større end 1e-12
		gsl_blas_ddot(u, grad_Fs, &uy); // prikprodukt u^T*y = uy
		gamma = uy/(2.0*sy); // gamma bestemmes 
	
		gsl_blas_daxpy(-gamma,Dx,u); // u=u-gamma*s
		
		gsl_blas_dger(1.0/sy,u,Dx,B); // B= B + u*s^T/s^T*y
		gsl_blas_dger(1.0/sy,Dx,u,B); // B= B + s*u^T/s^T*y
	}
	
	// opdaterer til det nye step
	gsl_vector_memcpy(x,y); 
    fx = fy;
    }
    printf("iterations = %d\n", iter);
	
    gsl_matrix_free(B);
    gsl_vector_free(grad_F);
    gsl_vector_free(grad_Fs);
    gsl_vector_free(y);  
    gsl_vector_free(By); 
    gsl_vector_free(Dx);
    gsl_vector_free(u);
    
}
		