#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include <stdio.h>

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
	gsl_matrix* QT = gsl_matrix_alloc(Q->size2,Q->size1);
	gsl_matrix_transpose_memcpy(QT, Q);
	gsl_blas_dgemv(CblasNoTrans, 1.0, QT, b, 0.0, x);
	
	int i=x->size-1; //size-1 to ensure proper indexing. i is the rows in R.
	for(i=i;i>=0;i--)
		{
		gsl_vector_set(x,i,(gsl_vector_get(x,i)-Rsum(R,x,i))/gsl_matrix_get(R,i,i));
		}
    gsl_matrix_free(QT);
}