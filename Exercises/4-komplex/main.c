#include"komplex.h"
#include"stdio.h"

int main(){
	komplex a = {2,4}, b = {5,8};

	printf("testing komplex_add\n");
	komplex r = komplex_add(a,b);
	komplex R = {7,12};
	komplex_print("a=", a);
	komplex_print("b=", b);
	komplex_print("a+b sould = ", R);
	komplex_print("a+b actually =", r);
	
	printf("testing komplex_sub\n");
	komplex p = komplex_sub(a,b);
	komplex P = {-3,-4};
	komplex_print("a=", a);
	komplex_print("b=", b);
	komplex_print("a-b sould = ", P);
	komplex_print("a-b actually =", p);
	
return 0;
}
