#include<stdio.h>
#include<limits.h>
#include<float.h>

int main(){
// opgave 1
//	int i=1; while(i+1>i) {
//		i++;
//		if(i > 1000000000) break;
//}
//	printf("my max int = %i\n", i);
//
//	i=1; for(i=1; i+1>i; i++) {
//		if(i > 1000000000) break;
//}
//	printf("my max int = %i\n", i);
//
//	i=1; do{printf("my max int = %i\n", i);
//		i++;
//		}while(i+1>1);
//

	double x=1; while(1+x!=1){x/=2;} x*=2;
	printf("x_double = %g\n", x);

	float y=1; while(1+y!=1){y/=2;} y*=2;
	printf("y_float = %f\n", y);

	long double z = 1; while(1+z!=1){z/=2;} z*=2;
	printf("z_longdouble = %Lg\n", z);

	double e; for(e=1; 1+e!=1; e/=2){} e*=2;
	printf("double with for loop = %g\n", e);

// opgave 2
	int max = INT_MAX/3;
	float sum_up_float=0.0;
	int i=1; for(i=1; i<=max; i++){sum_up_float += 1.0f/i;}

	float sum_down_float=0.0; 
	i=0; for(i=0; i<=max; i++){sum_down_float += 1.0f/(max-i);}

	printf("sum_up_float = %f\n", sum_up_float);
	printf("sum_down_float = %f\n", sum_down_float);

return 0;
}
