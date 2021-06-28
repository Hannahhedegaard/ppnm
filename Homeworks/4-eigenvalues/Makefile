CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: eigenfunctions.pyxplot.png timediagonalization.pyxplot.png out.txt

eigenfunctions.pyxplot.png: main eigenfunctions.txt Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Eigenfunctions";\
	plot \
	"eigenfunctions.txt" using 1:2 with line title "n=1"\
	,"eigenfunctions.txt" using 1:3 with line title "n=2"\
	,"eigenfunctions.txt" using 1:4 with line title "n=3"\
	,"exacteigenfunctions.txt" using 1:2 with line title "n=1 (exact)"\
	,"exacteigenfunctions.txt" using 1:3 with line title "n=2 (exact)"\
	,"exacteigenfunctions.txt" using 1:4 with line title "n=3 (exact)"\
	' | tee tmp.gpi | pyxplot

timediagonalization.pyxplot.png: main Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "N";\
	set ylabel "Time (s)";\
	set title "Task C";\
	a=0.00000058; b=3.1; f(x)=a*x**b ;\
	c=0.00000021; d=3.2; g(x)=c*x**d ;\
	plot \
	"out_time_jacobi.txt" using 1:2 with points title "Implementation"\
	,"out_time_gsl.txt" using 1:2 with points title "GSL"\
	,"out_time_jacobiop.txt" using 1:2 with points title "Optimized implementation"\
	,f(x) with lines title "fit: 5.8e-7*n**(3.1)"\
	,g(x) with lines title "fit: 2.1e-7*n**(3.2)"\
	' | tee tmp.gpi | pyxplot

eigenfunctions.txt: main

out.txt: main 
	./$< > $@

main: main.o functions.o

clean:
	$(RM) main main.o functions.o *.png *.txt *.gpi