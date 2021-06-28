#include<stdio.h>
#include<limits.h>
#include<float.h>

int main(){
	printf("Exercise 1. i)\n");
	printf("Finding maxiumum integer.\n");
	int i=0; while(i+1>i) {
		i++;
		}
	printf("Maximum integer with while loop = %i\n", i);

	i=0; for(i=1; i+1>i; i++) {
		}
	printf("Maximum integer with for loop = %i\n", i);

	i=0; do{i++;}
		while(i+1>1);
	printf("Maximum integer with do while loop = %i\n", i);

	printf("INT_MAX = %i\n",INT_MAX);
	
	printf("\nExercise 1. ii)\n");
	printf("Finding minimum integer.\n");
	
	i=0; while(i-1<i) {
		i++;
		}
	printf("Minimum integer with while loop = %i\n", i);

	i=0; for(i=1; i-1<i; i++) {
		}
	printf("Minimum integer with for loop = %i\n", i);

	i=0; do{
		i++;
		}while(i-1<i);
	printf("Minimum integer with do while loop = %i\n", i);

	printf("INT_MIN = %i\n",INT_MIN);
	
	printf("\nExercise 1. iii)\n");
	printf("Finding machine epsilon.\n");
	float y=1; while(1+y!=1){y/=2;} y*=2;
	printf("Machine eps for float with while loop = %.10g\n", y);
	y = 1; for(y=1; 1+y!=1; y/=2){} y*=2;
	printf("Machine eps for float with for loop = %.10g\n", y);
	y = 1; do{y/=2;}
	while(1+y!=1); y*=2;
	printf("Machine eps for float with do while loop = %.10g\n", y);
	printf("FLT_EPSILON = %.10g\n\n",FLT_EPSILON);
	
	double x=1; while(1+x!=1){x/=2;} x*=2;
	printf("Machine eps for double with while loop = %.10g\n", x);
	x = 1; for(x=1; 1+x!=1; x/=2){} x*=2;
	printf("Machine eps for double with for loop = %.10g\n", x);
	x = 1; do{x/=2;}
	while(1+x!=1); x*=2;
	printf("Machine eps for double with do while loop = %.10g\n", x);
	printf("DBL_EPSILON  = %.10g\n\n",DBL_EPSILON);
	
	long double z = 1; 
	while(1+z!=1){z/=2;} 
	z*=2;
	printf("Machine eps for long double with while loop = %.10Lg\n", z);
	long double g; 
	for(g=1; 1+g != 1; g/=2){} 
	g*=2;
	printf("Machine eps for long double with for loop = %.10Lg\n", g);
	long double c = 1; 
	do{c/=2;}
	while(1+c!=1); 
	c*=2;
	printf("Machine eps for long double with do while loop = %.10Lg\n", c);
	printf("LDBL_EPSILON = %.10Lg\n\n", LDBL_EPSILON);
	
	printf("\nExercise 2. i)\n");
	printf("Calculating sum.\n");
	int max = INT_MAX/3;
	float sum_up_float=0;
	i=1; for(i=1; i<=max; i++){
		sum_up_float += 1.0f/i;
		}

	float sum_down_float=0; 
	for(int i=max; i>0; i--){sum_down_float += 1.0f/i;}

	printf("sum_up_float = %f\n", sum_up_float);
	printf("sum_down_float = %f\n", sum_down_float);

	printf("\nExercise 2. ii)\nexplain the difference\n");
	printf("The 'up' sum is smaller than the 'down sum.\n");
	printf("This is because the 'up' sum has small terms at the end of the sum,\n");
	printf("which are outside of the machine precision and will therefore not count.\n");
	printf("In the 'down' sum, the small terms at the beginning will count.\n");
	
	printf("\nExercise 2. iii)\ndoes the sum converge?\n");
	printf("No, harmonic series are divergent\n");
	
	printf("\nExercise 2. iv)\nCalculate sums with double\n");
	double sum_up_double=0;
	i=1; for(i=1; i<=max; i++){
		sum_up_double += 1.0f/i;
		}

	double sum_down_double=0; 
	for(int i=max; i>0; i--){sum_down_double += 1.0f/i;}
	
	printf("sum_up_double = %f\n", sum_up_double);
	printf("sum_down_double = %f\n", sum_down_double);
	
	printf("The sums are now the same because double precision is smaller than float precision.");
	
return 0;
}
