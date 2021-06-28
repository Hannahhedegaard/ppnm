#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>

void plainmcA(double f(gsl_vector* x),gsl_vector* a,gsl_vector* b,int N,gsl_vector* results);
void plainmcB(double f(gsl_vector* x),gsl_vector* a,gsl_vector* b,int N,gsl_vector* results);

double f1(gsl_vector * x){
	return sqrt(gsl_vector_get(x,0));
	}

double f2(gsl_vector *x){
	if (1-cos(gsl_vector_get(x,0))*cos(gsl_vector_get(x,1))*cos(gsl_vector_get(x,2))==0) return 10;
	return 1/pow(M_PI,3)*pow(1-cos(gsl_vector_get(x,0))*cos(gsl_vector_get(x,1))*cos(gsl_vector_get(x,2)),-1);
	}

int main(){
	printf("Task A and B: Plain and quasi-random Monte Carlo integration.\n");
	printf("The two functions are implemented in functions.c.\n");
	int dim1=1;
	gsl_vector*a1=gsl_vector_alloc(dim1);
	gsl_vector*b1=gsl_vector_alloc(dim1);
	gsl_vector_set(a1,0,0);
	gsl_vector_set(b1,0,1);
	int N1=100000;
	gsl_vector* results1=gsl_vector_alloc(2);
    plainmcA(f1,a1,b1,N1,results1);
	printf("The functions are checked by calculating:\n");
    printf("Integral of sqrt(x) from 0 to 1 with plain MC        = %.10g, error = %.10g\n",gsl_vector_get(results1,0),gsl_vector_get(results1,1));
    plainmcB(f1,a1,b1,N1,results1);
    printf("Integral of sqrt(x) from 0 to 1 with quasi-random MC = %.10g, error = %.10g\n",gsl_vector_get(results1,0),gsl_vector_get(results1,1));
    printf("Exact result: %.10g\n", 2.0/3); 
    printf("Comparison of the errors for the two methods are shown in Figure 1.\n\n"); 
	FILE * file1=fopen("plain_int1.txt","w");
	FILE * file2=fopen("qrandom_int1.txt","w");
	for(int N=10000; N<100000; N+=5000){
		plainmcA(f1,a1,b1,N,results1);
		fprintf(file1,"%d %.25g %.25g\n",N,gsl_vector_get(results1,0),gsl_vector_get(results1,1));
		
		plainmcB(f1,a1,b1,N,results1);
		fprintf(file2,"%d %.25g %.25g\n",N,gsl_vector_get(results1,0),gsl_vector_get(results1,1));
		}
	fclose(file1);
	fclose(file2);
	
	printf("The functions are further tested by calculating the integral given in Task A.\n");

	int dim2=3;
	gsl_vector*a2=gsl_vector_alloc(dim2);
	gsl_vector*b2=gsl_vector_alloc(dim2);
	gsl_vector_set(a2,0,0);
	gsl_vector_set(b2,0,M_PI);
	gsl_vector_set(a2,1,0);
	gsl_vector_set(b2,1,M_PI);
	gsl_vector_set(a2,2,0);
	gsl_vector_set(b2,2,M_PI);
	int N2=100000;
	
	gsl_vector* results2=gsl_vector_alloc(2);
	plainmcA(f2,a2,b2,N2,results2);
	printf("Integral of integral in task A with plain MC        = %.10g, error = %.10g\n",gsl_vector_get(results2,0),gsl_vector_get(results2,1));
	plainmcB(f2,a2,b2,N2,results2);
    printf("Integral of integral in task A with quasi-random MC = %.10g, error = %.10g\n",gsl_vector_get(results2,0),gsl_vector_get(results2,1));
    printf("Exact result: 1.3932039296856768591842462603255\n"); 
	printf("Comparison of the errors for the two methods are shown in Figure 2.\n\n"); 
	
	FILE * file3=fopen("plain_int2.txt","w");
	FILE * file4=fopen("qrandom_int2.txt","w");
	for(int N=10000; N<100000; N+=2500){
		plainmcA(f2,a2,b2,N,results2);
		fprintf(file3,"%d %.25g %.25g\n",N,gsl_vector_get(results2,0),gsl_vector_get(results2,1));
		
		plainmcB(f2,a2,b2,N,results2);
		fprintf(file4,"%d %.25g %.25g\n",N,gsl_vector_get(results2,0),gsl_vector_get(results2,1));
		}
	fclose(file3);
	fclose(file4);
	
	printf("Comment: the two figures show, that the errors are smaller for the quasi-random method.\n");
	printf("The errors for the integral from task A with the plain method are widely distributed, and a fit is hard to make.\n");
	printf("Nevertheless, the errors are clearly smaller for the quasi-random method.\n");
	
return 0;
}
