#include <math.h>
#include <float.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>

static const double DELTA=1.0/524288;
void numeric_gradient
(double F(gsl_vector*), gsl_vector*x, gsl_vector*grad){
	double fx=F(x);
	for(int i=0;i<x->size;i++){
		double dx,xi=gsl_vector_get(x,i);
		if(fabs(xi)<sqrt(DELTA)) dx=DELTA;
		else dx=fabs(xi)*DELTA;
		gsl_vector_set(x,i,xi+dx);
		gsl_vector_set(grad,i,(F(x)-fx)/dx);
		gsl_vector_set(x,i,xi);
	}
}

void GS_decomp(gsl_matrix* A, gsl_matrix* R);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

void matrix_print(char s[], gsl_matrix* M);
void vector_print(char s[], gsl_vector* v);


void qnewton(double f(gsl_vector* x),gsl_vector* x, double eps){
	int n = x->size;
	int iter = 0;
	double lambda;
	double fx;
	double fy;
	double uy;
	double sy;
	double s_grad;
	const double dx = 1.0/524288;//sqrt(DBL_EPSILON);

	gsl_matrix* B = gsl_matrix_alloc(n,n);
	gsl_matrix* as = gsl_matrix_alloc(n,n);
	gsl_matrix* sa = gsl_matrix_alloc(n,n);
	gsl_matrix* dB = gsl_matrix_alloc(n,n);
	gsl_matrix_set_identity(B);
    gsl_vector* grad_F = gsl_vector_alloc(n);         
    gsl_vector* grad_Fs = gsl_vector_alloc(n);         
    gsl_vector* Dx = gsl_vector_alloc(n);    
    gsl_vector* y = gsl_vector_alloc(n);            
    gsl_vector* By = gsl_vector_alloc(n);    
    gsl_matrix* a = gsl_matrix_alloc(n,n);    
    gsl_matrix* s = gsl_matrix_alloc(n,n);    
    gsl_vector* u = gsl_vector_alloc(n);    
    gsl_matrix* R = gsl_matrix_alloc(n,n);
   
    while(1){
	iter += 1;
	if (iter>2000) break;
        fx = f(x); 

	numeric_gradient(f,x,grad_F);
	
	if(gsl_blas_dnrm2(grad_F)<eps){
		fprintf(stderr,"qnewton: |grad|<eps\n");
		printf("break1\n");
		break;
		}

	gsl_blas_dgemv(CblasNoTrans, -1, B, grad_F, 0, Dx);

    if(gsl_blas_dnrm2(Dx) < dx) {
		fprintf(stderr,"qnewton: |Dx|=%g < dx=%g\n",gsl_blas_dnrm2(Dx),dx);
		printf("break2\n");
        break;
        }

        lambda = 1;
	while(1){
		gsl_vector_memcpy(y,x);
		gsl_vector_add(y,Dx);
		fy = f(y);
		gsl_blas_ddot(Dx, grad_F, &s_grad);
		if(fy<fx+1e-4*s_grad){
			fprintf(stderr,"qnewton: fy=%g < fx=%g\n",fy,fx);
			break;
		}
		if (lambda < 1.0/32) {
			fprintf(stderr,"qnewton: lambda<dx !!!\n");
			gsl_matrix_set_identity(B);
			break;
		}
            lambda/=2;
            gsl_vector_scale(Dx,0.5);
        }

	// Opdatering af B til B + Î´B
	gsl_vector_add(Dx, x); 	// s <- s + x

	numeric_gradient(f,Dx,grad_Fs);
	
	// u og y-vektor bestemmes
	gsl_vector_sub(Dx, x); // Vektoren Dx er lig s. (s <- s - x)	
	gsl_vector_memcpy(u,Dx); // u = s
	gsl_vector_sub(grad_Fs, grad_F); // y-vektor bestemmes y = grad_Fs
	gsl_blas_dgemv(CblasNoTrans, 1.0, B, grad_Fs, 0.0, By); // B*y -> By
	gsl_vector_sub(u, By); // u = s - By
	// prikprodukter uy og sy bestemmes
	gsl_blas_ddot(u, grad_Fs, &uy); // prikprodukt u^T*y = uy
	gsl_blas_ddot(Dx, grad_Fs, &sy); // prikprodukt s^T*y = sy
	
	if (fabs(sy)>eps){ // laver B om hvis sy er stÃ¸rre end eps
		gsl_blas_dger(1/sy,u,Dx,B);
	}
	
	// opdaterer til det nye step
	gsl_vector_memcpy(x,y); 
    fx = fy;
    }
    printf("Number of steps = %d\n", iter);
	
    gsl_matrix_free(B);
    gsl_matrix_free(as);
    gsl_matrix_free(sa);
    gsl_matrix_free(dB);
    gsl_vector_free(grad_F);
    gsl_vector_free(grad_Fs);
    gsl_vector_free(y);  
    gsl_vector_free(By); 
    gsl_vector_free(Dx);
    gsl_matrix_free(a);
    gsl_matrix_free(s);
    gsl_vector_free(u);
    gsl_matrix_free(R);
}

		