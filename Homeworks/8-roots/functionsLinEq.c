#include <math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_linalg.h>

void matrix_print(char s[], gsl_matrix* M){
    printf("%s\n",s);
	for(int i=0;i< M->size1;i++){
        for(int j = 0; j< M->size2; j++){
            //printf("i, j = %10d %10d \n", i,j);
        printf("%10g ", gsl_matrix_get(M,i,j));}
    printf("\n");}
}

void GenerateMatrix(gsl_matrix* A, int N, int M) {
    
    for(int i = 0; i<N; i++) {
        for(int j = 0; j<M; j++) {
            gsl_matrix_set(A, i, j, (double)rand()/RAND_MAX);}   
    }
    
}

void GS_decomp(gsl_matrix* A, gsl_matrix* R) {
   int N = A->size1; //Size1 er antal rækker
   int M = A->size2;   //Size1 er antal søjler
   for (int i =0; i<M; i++){
       
       //Sets Rii=sqrt(ai^T*ai)
       gsl_vector* ai = gsl_vector_alloc(N); 
       gsl_matrix_get_col(ai, A, i);
       double dot; 
       gsl_blas_ddot(ai, ai, &dot); 
       gsl_matrix_set(R, i, i, sqrt(dot));
       
       //Sets qi = ai/Rii
       gsl_vector_scale(ai, (double)1/sqrt(dot)); 
       gsl_matrix_set_col(A, i, ai);
       
       
       for (int j = i+1; j < M; j++){ 
            //Sets Rij = qi^T*aj
            double dot2; 
            gsl_vector* aj = gsl_vector_alloc(N); 
            gsl_matrix_get_col(aj, A, j);
            //gsl_matrix_get_col(ai, A, i);
            gsl_blas_ddot(ai, aj, &dot2); 
            gsl_matrix_set(R, i, j, dot2);
            
            //Sets aj = aj-qi*Rij
            gsl_vector_scale(ai,dot2);
            gsl_vector_sub(aj, ai); 
            gsl_matrix_set_col(A, j, aj); //aj = aj-qi*Rij
            gsl_vector_scale(ai,(double)1/dot2);
            gsl_vector_free(aj);
            
            }
    gsl_vector_free(ai);    
   }
}

double Rsum(gsl_matrix *R, gsl_vector *x, int i){
	if(i==x->size-1) return 0;
    
	double sum=0;
	int j=x->size-1; //-1 to ensure proper indexing
	for(j=j;j>i;j--)
		{
		sum+=gsl_matrix_get(R,i,j)*gsl_vector_get(x,j);
		}
return sum;
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){

	gsl_matrix_transpose(Q);
	gsl_blas_dgemv(CblasNoTrans, 1.0, Q, b, 0.0, x);

	int i=x->size-1; //size-1 to ensure proper indexing. i is the rows in R.
	for(i=i;i>=0;i--)
		{
		gsl_vector_set(x,i,(gsl_vector_get(x,i)-Rsum(R,x,i))/gsl_matrix_get(R,i,i));
		}
    gsl_matrix_transpose(Q);
}
 
void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
	gsl_matrix_set_identity(B); 
	for(int i=0; i<Q->size2; i++){
		gsl_vector_view v = gsl_matrix_column(B,i);
        gsl_vector* v_save = gsl_vector_alloc(B->size1);
        gsl_vector_memcpy(v_save,&v.vector);
		GS_solve(Q, R, v_save, &v.vector);
        gsl_vector_free(v_save);
	}	
}