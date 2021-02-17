#include<stdio.h>
#include<math.h>
int main(int argc, char** argv){
	double x;
	int items;
	FILE* my_in_stream = fopen(argv[1], "r");
	FILE* my_out_stream = fopen(argv[2], "w");
	
	do{
		items=fscanf(my_in_stream, "%lg", &x);
		fprintf(my_out_stream, "x=%g cos(x)=%g\n", x, cos(x));
	}while(items!=EOF);
fclose(my_in_stream);
fclose(my_out_stream);
return 0;
}
