#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++)printf("%10g \n",gsl_vector_get(v,i));
//	printf("\n");
}



void rkstep12(void f(double x,gsl_vector* y,gsl_vector* dydx), double x, gsl_vector* yx,        
	double h, gsl_vector* yh, gsl_vector* dy){
    int i; 
    int n = yx->size;
    gsl_vector* k0 = gsl_vector_alloc(n); 
    gsl_vector* yt = gsl_vector_alloc(n);
    gsl_vector* k12 = gsl_vector_alloc(n);
    f(x, yx, k0);     
    for(i=0; i<n; i++){ 
        gsl_vector_set(yt,i,gsl_vector_get(yx,i)+gsl_vector_get(k0,i)*h/2);}
    f(x+h/2,yt,k12); 
    for(i=0; i<n; i++){
        gsl_vector_set(yh,i,gsl_vector_get(yx,i)+gsl_vector_get(k12,i)*h);}
    for(i=0; i<n; i++){
        gsl_vector_set(dy,i,(gsl_vector_get(k0,i)-gsl_vector_get(k12,i))*h/2);
    }    
}
    
void driver(void f(double x,gsl_vector* y ,gsl_vector* dydx), double a, gsl_vector* ya,               
	double b, gsl_vector* yb, double h, double acc, double eps, char filename[]){ 
	FILE * file = fopen(filename, "w"); 
	int n = ya->size;
	gsl_vector* dy = gsl_vector_alloc(n); 
	double err = 0, normy = 0, tol = 0; 
	int i; 
	while(a < b){
		if(a+h>b) {
			h = b-a; }	
		rkstep12(f, a, ya, h, yb, dy); 
		double s = 0.0; 
		for (i=0; i<n; i++) {
			s+= gsl_vector_get(dy,i)*gsl_vector_get(dy, i); 
			err = sqrt(s); 
		}
		s = 0.0; 
		for (i=0; i<n; i++) {
			s+= gsl_vector_get(yb,i)*gsl_vector_get(yb, i); 
			normy = sqrt(s); 
		}
    
		tol = (normy*eps+acc)*sqrt(h/(b-a)); 
		if(err<tol){
			a += h; 
			for(i = 0; i<n; i++) {
				gsl_vector_set(ya,i, gsl_vector_get(yb,i)); 
			}
            fprintf(file, "%10g ", a); 
			for (i = 0; i<ya->size; i++) {
				fprintf(file, "%10g ", gsl_vector_get(ya,i)); 
			}
        fprintf(file, "\n"); 
		}
    
		if(err>0){
			h*=pow(tol/err, 0.25)*0.95; 
			}
		else{h*=2;}
	}
fclose(file);
}






