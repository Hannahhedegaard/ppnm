#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_blas.h>

void vector_print_2(char s1[], char s2[], gsl_vector* v1, gsl_vector* v2){
	printf("%s    %s\n",s1, s2);
	for(int i=0;i<v1->size;i++)printf("%8g %12g\n",gsl_vector_get(v1,i), gsl_vector_get(v2,i));
}

typedef struct {gsl_vector *x, *y , *b, *c, *d;} akima_spline;
akima_spline* akima_spline_alloc(gsl_vector *x, gsl_vector *y);
double akima_spline_eval(akima_spline *s, double z);
void akima_spline_free(akima_spline *s);

typedef struct {gsl_vector *x, *y, *b, *c, *d;} cubic_spline;
cubic_spline* cubic_spline_alloc(gsl_vector *x, gsl_vector* y);
double cubic_spline_eval(cubic_spline *s, double z);
void cubic_spline_free(cubic_spline *s);

int main(){
	printf("Akima sub-spline.\n");
	printf("The Akima sub-spline is implemented in functions.c.\n");
	printf("The implementation is now tested on two different set of data-points.\n");
	printf("---------------------------------------------------\n");
	printf("Dataset 1: see input1.dat for original data-points.\n");
	printf("This will lead to a replication of Figure 1.2 from the book.\n\n");
	
	//The length of the dataset is saved in the integer N1.
	int N1 = 0; int ch;
	FILE* fp = fopen("input1.dat","r");
	while(!feof(fp)){ch = fgetc(fp); if(ch == '\n'){N1++;}}
	fclose(fp);
	printf("Length of dataset 1: %d\n", N1);
	
	//Memory is allocated to two vectors x and y, which will contain the x and y values for each point.
	gsl_vector * x1 = gsl_vector_alloc(N1);	
	gsl_vector * y1 = gsl_vector_alloc(N1);

	//The data-points are loaded in to the two vectors x and y.
	double i; double j;
	int items; int n=0;
	FILE* input1 = fopen("input1.dat","r");
	while((items=fscanf(input1,"%lg %lg",&i,&j))!=EOF){
		gsl_vector_set(x1, n, i);
		gsl_vector_set(y1, n, j);
		n++;
		}	
	vector_print_2("x-values:", "y-values:", x1, y1);
	fclose(input1);
	
	FILE * f1 = fopen("data1.txt","w"); //data-file: contains the spline data.
	double z1 = gsl_vector_get(x1,0); //z (used in the eval-functions) are determined as the first x-point.
	akima_spline *s1 = akima_spline_alloc(x1,y1); //the Akima sub-spline is allocated.
	cubic_spline *sc1 = cubic_spline_alloc(x1,y1); //the cubic spline is allocated - for comparison with the akima sub-spline.
	//The two splines are calculated for the points z, and written into the data1.txt file.
	for((z1=z1); z1<=gsl_vector_get(x1,N1-1); z1+=1./16){
		double eval = akima_spline_eval(s1,z1);
		double eval_cubic = cubic_spline_eval(sc1,z1);
		fprintf(f1,"%g %g %g\n",z1,eval, eval_cubic);
		}
	fclose(f1);
	//Free the allocated memory:
	akima_spline_free(s1);
	cubic_spline_free(sc1);
	
	printf("\nThe Akima sub-spline is calculated, see results in data1.txt.\n");
	printf("Furthermore, the cubic spline is calculated, results also in data1.txt.\n");
	printf("The results are shown in Figure1.pyxplot.png, and as this is a replication of Figure 1.2 in the book,\n");
	printf("the two figures can be compared.\n");
	printf("\nComment: the two figures are very similar - and it is clear that the Akima sub-spline\n");
	printf("wiggels less than the cubic spline, as expected.\n");
	
	printf("---------------------------------------------------\n");
	printf("Dataset 2: see input2.dat for original data-points.\n");
	printf("This data-set is made to show how the Akima sub-spline interpolates,\n");
	printf("when one of the data-points lies away from the others.\n\n");
	
	//The length of the dataset is saved in the integer N2.
	int N2=0; int ch2;
	FILE* fp2 = fopen("input2.dat","r");
	while(!feof(fp2)){ch2 = fgetc(fp2); if(ch2 == '\n'){N2++;}}
	fclose(fp2);
	printf("Length of dataset 2: %d\n", N2);
	
	//Memory is allocated to two vectors x and y, which will contain the x and y values for each point.
	gsl_vector * x2 = gsl_vector_alloc(N2);	
	gsl_vector * y2 = gsl_vector_alloc(N2);

	//The data-points are loaded in to the two vectors x and y.
	double k; double l;
	int items2; int m=0;
	FILE* input2 = fopen("input2.dat","r");
	while((items2=fscanf(input2,"%lg %lg",&k,&l))!=EOF){
		gsl_vector_set(x2, m, k);
		gsl_vector_set(y2, m, l);
		m++;
		}
	vector_print_2("x-values:", "y-values:", x2, y2);
	fclose(input2);
	
	FILE * f2 = fopen("data2.txt","w");
	double z2 = gsl_vector_get(x2,0);
	akima_spline *s2 = akima_spline_alloc(x2,y2); //again is both the Akima sub-spline and cubic spline calculated.
	cubic_spline *sc2 = cubic_spline_alloc(x2,y2);
	//The two splines are calculated for the points z, and written into the data2.txt file.
	for((z2=z2); z2<=gsl_vector_get(x2,N2-1); z2+=1./16){
		double eval2 = akima_spline_eval(s2,z2);
		double eval_cubic2 = cubic_spline_eval(sc2,z2);
		fprintf(f2,"%g %g %g\n",z2,eval2, eval_cubic2);
		}
	fclose(f2);
	//Free the allocated memory:
	akima_spline_free(s2);
	cubic_spline_free(sc2);
	
	printf("\nThe Akima sub-spline is calculated, see results in data2.txt.\n");
	printf("Furthermore, the cubic spline is calculated, results also in data2.txt.\n");
	printf("The results are shown in Figure2.pyxplot.png.\n");
	
	printf("\nComment: It is shown in the figure, that the Akima sub-spline makes a smooth interpolation\n");
	printf("through the data-points, even though there is an outlier. The cubic spline, on the other hand,\n");
	printf("shows the typical wiggels.\n");
	
	return 0;
}