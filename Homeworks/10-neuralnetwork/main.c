#include <math.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "ann.h"

void qnewton(double f(gsl_vector* x),gsl_vector* x,double eps);

void vector_print(char s[], gsl_vector* v){
	printf("%s",s);
	for(int i=0;i<v->size;i++)printf("%10g \n",gsl_vector_get(v,i));
	printf("\n");
}
	
double activationfunction(double x){
	return x*exp(-(x*x));
}

double fittingfunction(double x){
	return cos(5*x-1)*exp(-fabs(x));
}

gsl_vector * xfunc;
gsl_vector * yexact;
ann * network;

double cost_function(gsl_vector* p){
		gsl_vector_memcpy(network->params,p);
		double sum=0;
		for(int i=0;i<xfunc->size;i++){
			double xi=gsl_vector_get(xfunc,i);
			double yi=gsl_vector_get(yexact,i);
			double fi=ann_response(network,xi);
			sum+=fabs(fi-yi);
		}
		return sum/xfunc->size;
}

int main(){
	printf("Task A\n");
	printf("In this exercise we construct a simplest artificial neural network which will be trained to interpolate a tabulated function.\n");
	printf("The activationfunction is a Gaussian wavelet, f(x) = x*exp(-x^2).\n");
	int n=5; //number of neurons
	network=ann_alloc(n,*activationfunction);

	double a=-3;
	double b=3;
	
	for(int i=0;i<n;i++){
		gsl_vector_set(network->params,i*3,a+(b-a)*i/(network->n-1));
		gsl_vector_set(network->params,i*3+1,1);
		gsl_vector_set(network->params,i*3+2,1);
	}
	vector_print("network->params\n",network->params);
	int np=100;
	xfunc=gsl_vector_alloc(np);
	yexact=gsl_vector_alloc(np);
	for (int i=0;i<np;i++){
		double x=a+(b-a)*i/(np-1);
		double f=fittingfunction(x);
		gsl_vector_set(xfunc,i,x);
		gsl_vector_set(yexact,i,f);
	}
	
	printf("The fitting function is cos(5*x-1)*exp(-|x|)\n"); 
	printf("Training...\n");
	ann_train(network,xfunc,yexact);
	printf("The final parameters:\n"); 
	vector_print("network->params\n",network->params);
	printf("A comparison of the estimated and exact y-values is shown in Figure 1.\n");
	
	gsl_vector * ygaet=gsl_vector_alloc(np);
	
	for(int i=0;i<xfunc->size;i++){
		gsl_vector_set(ygaet,i,ann_response(network,gsl_vector_get(xfunc,i)));
	}

	FILE* f=fopen("data.txt","w");
	for(int i=0;i<np;i++){
		fprintf(f,"%g %g %g\n",gsl_vector_get(xfunc,i),gsl_vector_get(yexact,i),gsl_vector_get(ygaet,i));
	}
	printf("-------------------------------\n");
	printf("Task B\n");
    printf("The method is now modified so the network can also approximate derivative and anti-derivative.\n");
    printf("Figure 2 shows the comparison of the estimated and exact y-values, together with the derivative and anti-derivative.\n");
	double deltax=1.0/16;
	FILE* f2=fopen("dataB.txt","w");
	for(double step=a;step<=b;step+=deltax){
		double gradient=ann_gradient(network,step);
		double value=ann_response(network,step);
		double integral=ann_integral(network,step);
		double exactvalue=fittingfunction(step);
		fprintf(f2,"%g %g %g %g %g\n",step,gradient,value,integral,exactvalue);
	}
	
return 0;
}