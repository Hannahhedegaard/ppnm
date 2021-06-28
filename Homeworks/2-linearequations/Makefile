CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

Decomp_time.png: out_decomp.txt out_gsl.txt Makefile main out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "N (Matrix of size NxN)";\
	set ylabel "Time (s)";\
	set title "Performance of the function compared to GSL: QR-Decomposition";\
    a=0.0000000167; b=3.0; f(x)=a*x**b ;\
	c=0.0000000066; d=3.0; g(x)=c*x**d ;\
	plot \
	"out_decomp.txt" using 1:2 with points title "Implemented decomposition"\
	,"out_gsl.txt" using 1:2 with points title "GSL decomposition"\
    ,f(x) with lines title "fit 10.7e-9*n**(3.0)"\
	,g(x) with lines title "fit 5.3e-9*n**(3.0)"\
	' | pyxplot

out.txt: main 
	./$< > $@

out_decomp.txt out_gsl.txt: main
        
main: functions.o main.o 

clean: 
	$(RM) main *.o *.txt *.gpi *.png *.log
