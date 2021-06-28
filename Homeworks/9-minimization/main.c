#include <math.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++)printf("%10g \n",gsl_vector_get(v,i));
}

void qnewton(double f(gsl_vector* x),gsl_vector* x, double eps);

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
	
	double mass = gsl_vector_get(x,0);
	double gamma = gsl_vector_get(x,1);
	double A = gsl_vector_get(x,2);
	
	double sum = 0;
	for (int i = 0; i<30; i++){
		sum += pow(A/(pow((gsl_vector_get(E,i)-mass),2)+pow(gamma,2)/4.0)-gsl_vector_get(sigma,i),2)/pow(gsl_vector_get(dsigma,i),2);
	}
return sum;
}

int downhill_simplex(double F(double *), double** simplex, int d, double simplex_size_goal);

int main(){
	printf("Task A: Quasi-Newton method with numerical gradient, back-tracking linesearch, and rank-1 update.\n");
	printf("The quasi-Newton minimization method is implemented in functions.c.\n");
	printf("The routine is checked by finding minimum in four different functions:\n");
	
	printf("\nFind minimum of function: f(x) = (x+1)^2 - 3.\n");
	int n = 1;
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector_set(x,0,-2);
	double eps = 1e-3;
	qnewton(f1, x, eps);
	vector_print("Calculated minimum =", x);
	printf("Exact minimum = -1\n");
	
	printf("\nFind minimum of function: f(x,y) = x^2 + y^2.\n");
	n = 2;
	gsl_vector* x1 = gsl_vector_alloc(n);
	gsl_vector_set(x1,0,-2);
	gsl_vector_set(x1,1,-2);
	qnewton(f2, x1, eps);
	vector_print("Calculated minimum =", x1);
	printf("Exact minimum = (0,0).\n");
	
	printf("\nFind minimum of Rosenbrock's valley function: f(x,y) = (1-x)^2 + 100(y-x^2)^2.\n");
	n = 2;
	gsl_vector* x2 = gsl_vector_alloc(n);
	gsl_vector_set(x2,0,18);
	gsl_vector_set(x2,1,15);
	qnewton(fRose, x2, eps);
	vector_print("Calculated minimum =", x2);
	printf("Exact minimum = (1,1).\n");
	
	printf("\nFind minimum of Himmelblau's function: f(x,y) = (x^2+y-11)^2+(x+y^2-7)^2.\n");
	n = 2;
	gsl_vector* x3 = gsl_vector_alloc(n);
	gsl_vector_set(x3,0,10);
	gsl_vector_set(x3,1,10);
	qnewton(fHimmel, x3, eps);
	vector_print("Calculated minimum =", x3);
	printf("Exact minimum = (3,2), (-2.805118, 3.131312), (-3.779310, -3.283186) or (3.584428, -1.848126).\n");
	
	printf("-----------------------------------------\n");
	printf("Task B: Higgs discovery.\n");
	printf("The Breit-Wigner function is fitted to the data, and the mass, width, and scale factor is determined.\n");
	n = 3;
	gsl_vector* x4 = gsl_vector_alloc(n);
	gsl_vector_set(x4,0,120);
	gsl_vector_set(x4,1,2);
	gsl_vector_set(x4,2,10);
	qnewton(fBW, x4, eps);
	vector_print("Calculated mass, width, and scale factor: =", x4);
	
	printf("-----------------------------------------\n");
	printf("Task C: Implement the downhill simplex method.\n");
	printf("The method is implemented in functions.c.\n");
	printf("The method is checked by calculating the minimum of the Rosenbrock function:\n");
	int d = 2;
	double* simplex[] = {
      (double[]){0.0, 0.0},
      (double[]){0.0, 5.0},
      (double[]){5.0, 0.0}
	};
	double sizegoal = 1e-5;
	int k; //antal steps
	k = downhill_simplex(f, simplex, d, sizegoal);
	printf("Number of steps = %d\n", k);
	
	printf("Compared to the method in task A, the downhill simplex method is more accurate and efficient.\n");
	
	return 0;
}