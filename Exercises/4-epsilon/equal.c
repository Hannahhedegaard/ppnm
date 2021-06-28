#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int equal(double a, double b, double tau, double epsilon){
	double i = fabs(a-b);
	if(i < tau){return 1;}

	if(i/(fabs(a)+fabs(b)) < epsilon/2){return 1;}

return 0;
}

int main(){
	printf("\n\nExercise 3.\nequal function\n");
	int x = equal(2.5, 2.6, 0.01, 0.01);
	printf("%d\n", x);
	int y = equal(2.5, 2.6, 0.11, 0.11);
	printf("%d\n", y);
return 0;
}
