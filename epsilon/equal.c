#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int equal(double a, double b, double tau, double epsilon){
	double i = fabs(a-b);
	if(i < tau){printf("1\n");}

	if(i/(fabs(a)+fabs(b)) < epsilon/2){printf("1\n");}

return 0;
}

int main(){
	equal(2.5, 2.6, 0.2, 0.2);
return 0;
}
