#include <math.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++)printf("%10g \n",gsl_vector_get(v,i));
//	printf("\n");
}

void qnewton(
	double f(gsl_vector* x), /* objective function */
	gsl_vector* x, /* on input: starting point, on exit: approximation to root */
	double eps /* accuracy goal, on exit |gradient| should be <eps */
	);

double f1(gsl_vector* x){
	double xx = gsl_vector_get(x,0);
    return pow((xx+1),2)-3;
}

double f2(gsl_vector* x){
	double xx = gsl_vector_get(x,0);
	double yy = gsl_vector_get(x,1);
	return (pow(xx,2)+pow(yy,2));
}

double fRose(gsl_vector* x){
	double xx = gsl_vector_get(x,0);
	double yy = gsl_vector_get(x,1);
	return (pow(1-xx,2)+100*pow((yy-pow(xx,2)),2));
}

double fHimmel(gsl_vector* x){
	double xx = gsl_vector_get(x,0);
	double yy = gsl_vector_get(x,1);
	return pow((xx*xx+yy-11),2)+pow((xx+yy*yy-7),2);
}

double f(double* x) {
   return pow(1 - x[0], 2) + 100*pow(x[1] - x[0]*x[0], 2);
}


gsl_vector* E;
gsl_vector* sigma;
gsl_vector* dsigma;

double fBW(gsl_vector* x){
	int N = 30;
	gsl_vector* E = gsl_vector_alloc(N);
	gsl_vector* sigma = gsl_vector_alloc(N);
	gsl_vector* dsigma = gsl_vector_alloc(N);

	FILE* file=fopen("opgB_data.txt","r");
	int m = 0;
	double i;
    double j;
	double k;
	int items;
    while((items=fscanf(file,"%lg %lg %lg",&i,&j,&k))!=EOF){
        gsl_vector_set (E, m, i);
        gsl_vector_set (sigma, m, j);
		gsl_vector_set (dsigma, m, k);
        m++;
        }
    fclose(file);
	
	//vector_print("E =", E);
	//vector_print("sigma =", sigma);
	//vector_print("dsigma =", dsigma);
	
	double mass = gsl_vector_get(x,0);
	double gamma = gsl_vector_get(x,1);
	double A = gsl_vector_get(x,2);
	
	//printf("mass = %g\n", mass);
	//printf("gamma = %g\n", gamma);
	//printf("A = %g\n", A);

	double sum = 0;
	for (int i = 0; i<30; i++){
		sum += pow(A/(pow((gsl_vector_get(E,i)-mass),2)+pow(gamma,2)/4.0)-gsl_vector_get(sigma,i),2)/pow(gsl_vector_get(dsigma,i),2);
	}
	//printf("sum = %g\n", sum);
	return sum;
}

int downhill_simplex(double F(double *), double** simplex, int d, double simplex_size_goal);

int main(){
	// opgave A
	int n = 1;
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector_set(x,0,-2);
	double eps = 1e-3;
	qnewton(f1, x, eps);
	vector_print("x =", x);
	
	n = 2;
	gsl_vector* x1 = gsl_vector_alloc(n);
	gsl_vector_set(x1,0,-2);
	gsl_vector_set(x1,1,-2);
	qnewton(f2, x1, eps);
	vector_print("x_f2 =", x1);
	
	n = 2;
	gsl_vector* x2 = gsl_vector_alloc(n);
	gsl_vector_set(x2,0,18);
	gsl_vector_set(x2,1,15);
	qnewton(fRose, x2, eps);
	vector_print("x_rose =", x2);
	
	n = 2;
	gsl_vector* x3 = gsl_vector_alloc(n);
	gsl_vector_set(x3,0,10);
	gsl_vector_set(x3,1,10);
	qnewton(fHimmel, x3, eps);
	vector_print("x_Himmel =", x3);
	
	// opgave B
	n = 3;
	gsl_vector* x4 = gsl_vector_alloc(n);
	gsl_vector_set(x4,0,120);
	gsl_vector_set(x4,1,2);
	gsl_vector_set(x4,2,10);
	qnewton(fBW, x4, eps);
	vector_print("x_BW =", x4);
	
	// Opgave C
	int d = 2;
	double* simplex[] = {
      (double[]){0.0, 0.0},
      (double[]){0.0, 5.0},
      (double[]){5.0, 0.0}
	};
	double sizegoal = 1e-5;
	//double fxs[d + 1];
	//double centroid[d];
	//int hi = 0;
	//int lo = 0;
	int k; //antal iterations
	k = downhill_simplex(f, simplex, d, sizegoal);
	printf("Iterations = %d\n", k);
	
	
	return 0;
}