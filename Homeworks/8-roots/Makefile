CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: Hydrogen.pyxplot.png C_firstcondition.pyxplot.png C_secondcondition.pyxplot.png C1andC2.pyxplot.png

Hydrogen.pyxplot.png: hydrogenatom.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Hydrogenatom";\
	plot \
	"hydrogenatom.txt" using 1:2 with points title "Calculated"\
    , x*exp(-x) with lines title "Exact Result"\
	' | tee tmp.gpi | pyxplot
    
    
C_firstcondition.pyxplot.png: hydrogenatom.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "rmax";\
	set ylabel "Energy";\
	set title "C1";\
	plot \
	"C_firstcondition.txt" using 1:2 with points title "Calculated"\
    , -0.5 with lines title "Exact Result"\
	' | tee tmp.gpi | pyxplot

    
C_secondcondition.pyxplot.png: hydrogenatom.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "rmax";\
	set ylabel "Energy";\
	set title "C2";\
	plot \
	"C_secondcondition.txt" using 1:2 with points title "Calculated"\
    , -0.5 with lines title "Exact Result"\
	' | tee tmp.gpi | pyxplot

C1andC2.pyxplot.png: hydrogenatom.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "rmax";\
	set ylabel "Energy";\
	set title "Comparisson of C1 and C2";\
	plot \
	"C_firstcondition.txt" using 1:2 with points title "C1"\
	,"C_secondcondition.txt" using 1:2 with points title "C2"\
    , -0.5 with lines title "Exact Result"\
	' | tee tmp.gpi | pyxplot



out.txt: main 
	./$< > $@
    
hydrogenatom.txt: main    

    
main: ODEfunctions.c functions.c main.c  functionsLinEq.c
   
    
clean: 
	$(RM) main *.o *.txt *.png *.gpi
    
    