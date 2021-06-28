CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: Figure1.pyxplot.png Figure2.pyxplot.png

Figure1.pyxplot.png: main Makefile plain_int1.txt qrandom_int1.txt out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "N";\
	set ylabel "error";\
	set title "error comparision";\
	a=0.235; b=-0.5; f(x)=a*x**b ;\
	c=1.3; d=-0.9; g(x)=c*x**d ;\
	plot \
	"plain_int1.txt" using 1:3 with points title "plain-method error "\
	,"qrandom_int1.txt" using 1:3 with points title "quasi-random-method error"\
	,f(x) with lines title "fit: 0.235*n**(-0.5)"\
	,g(x) with lines title "fit: 1.3*n**(-0.9)"\
	' | pyxplot
	
Figure2.pyxplot.png: main Makefile plain_int2.txt qrandom_int2.txt out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "N";\
	set ylabel "error";\
	set title "error comparision";\
	e=0.2; f=-0.2; i(x)=e*x**f ;\
	g=7.6; h=-0.9; j(x)=g*x**h ;\
	plot \
	"plain_int2.txt" using 1:3 with points title "plain-method error "\
	,"qrandom_int2.txt" using 1:3 with points title "quasi-random-method error"\
	, i(x) with lines title "fit: 0.2*n**(-0.2)"\
	, j(x) with lines title "fit: 7.6*n**(-0.9)"\
	' | pyxplot

plain_int1.txt:main

qrandom_int1.txt:main

plain_int2.txt:main

qrandom_int2.txt:main

out.txt: main 
	./$< > $@
    
main: functions.c main.c 
    
clean: 
	$(RM) main main.o functions.o *.txt *.png *.gpi