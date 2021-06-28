#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_blas.h>


void GenerateMatrix(gsl_matrix* A, int N, int M);
void matrix_print(char s[], gsl_matrix* M);
void GS_decomp(gsl_matrix* A, gsl_matrix* R);
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void vector_print(char s[], gsl_vector* v); 
void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);
void leastsq(gsl_vector* x, gsl_vector* y, gsl_vector* dy,gsl_vector* c,gsl_vector* dc, double f(int k, double z));
double f(int k, double z);

int main(){
    gsl_vector* t = gsl_vector_alloc(9);
    gsl_vector* y = gsl_vector_alloc(9);
    gsl_vector* dy = gsl_vector_alloc(9);
	gsl_vector* c = gsl_vector_alloc(2);
	gsl_vector* dc = gsl_vector_alloc(c->size);
	
	double i;
    double j;
	double k;
	double l;
	double m;
    int items;
    int n=0;
    FILE* my_in_stream=fopen("data_input.txt","r");
    while((items=fscanf(my_in_stream,"%lg %lg %lg %lg %lg",&i,&j,&k,&l,&m))!=EOF){
        gsl_vector_set (t, n, i);
        gsl_vector_set (y, n, j);
		gsl_vector_set (dy, n, m);
        n++;
        }
    fclose(my_in_stream);
	
	gsl_vector* y_ln = gsl_vector_alloc(9);
	gsl_vector_memcpy(y_ln,y);
	for (int i=0; i<y_ln->size; i++){
		gsl_vector_set(y_ln,i,log(gsl_vector_get(y_ln,i)));
	}
	
	leastsq(t, y_ln, dy, c, dc, f);
	
	FILE * f1=fopen("data_output.txt","w"); 
	for (double x = 1.0/16; x<15; x+=1.0/8) {
		double sum = 0;
		double upsum = 0;
		double downsum = 0;
		for (int k = 0; k<c->size; k++) {
			sum += gsl_vector_get(c,k)*f(k,x);
			upsum += (gsl_vector_get(c,k)+gsl_vector_get(dc,k))*f(k,x);
			downsum += (gsl_vector_get(c,k)-gsl_vector_get(dc,k))*f(k,x);
		}
		fprintf(f1, "%g %g %g %g\n", x, sum, upsum, downsum);
	}
	printf("Task A\nThe experimental data are plotted with errorbars and the best fit.\n");
	printf("The half-life time is calculated as: log(2)/λ, where λ is the slope from the fit.\n");
	printf("Half-life time = %g days\n", log(2)/gsl_vector_get(c,1));
	printf("Half-life time (today) = 3.6319 days\n"); //https://en.wikipedia.org/wiki/Isotopes_of_radium
	printf("\nTask B\nCalculating the uncertainties.\n");
	printf("The uncertainty is calculated from propagation of error: log(2)/λ^2 * σ_λ\n");
	printf("Uncertainty in half-life time = %g days\n", log(2)/pow(gsl_vector_get(c,1),2)*gsl_vector_get(dc,1)); //fejlophobningslov
	printf("Half-life time uncertainty (today) = 0.0023 days\n");
	printf("The modern uncertainty is a factor of 10 smaller than the one obtained in 1902.\nThis is expected as modern experiments are much more accurate.\n");
	printf("\nTask C\nTwo extra fits are plottet, where the fit coefficients are changed by δc.");
	
	gsl_vector_free(t);
    gsl_vector_free(y);
    gsl_vector_free(dy);
    gsl_vector_free(c);
    gsl_vector_free(y_ln);
return 0; 
}