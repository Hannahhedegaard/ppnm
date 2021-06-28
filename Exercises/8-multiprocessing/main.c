#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <stdlib.h>

typedef struct{
	double x;
	double y;}
	point;

typedef struct{int N; int Nin; unsigned int seed;} counts;

int in_out(point p){
	// returns 1 if point is in circle
	double x = p.x;
	double y = p.y;
	if(sqrt(pow(x,2)+pow(y,2))<1){return 1;}
		
return 0;
}

void* inside(void* arg){
	counts * param = (counts*) arg;
	int Nn=(*param).N;
	unsigned int seed=(*param).seed;
	double count = 0;
	for (int i = 0; i < Nn; i++){
	point p1 = {.x=((double)rand_r(&seed)/RAND_MAX), .y=((double)rand_r(&seed)/RAND_MAX)};
	if(in_out(p1)==1) count ++;
	else continue;
	}
	(*param).Nin=count;
return NULL;
}

int main(){
	int n = 1e6;
	counts nt1={.N=n/2, .Nin=0, .seed=1};
	counts nt2={.N=n/2, .Nin=0, .seed=2};
	pthread_t t1, t2;
	pthread_attr_t* attributes = NULL;
	pthread_create( &t1, attributes, inside, (void*)&nt1);
	pthread_create( &t2, attributes, inside, (void*)&nt2);
	pthread_join(t1, NULL);
	pthread_join(t2, NULL);
	double Nin=nt1.Nin+nt2.Nin;
	printf("π-approks = %.10g\n", 4*((double)Nin/n));
return 0;
}
