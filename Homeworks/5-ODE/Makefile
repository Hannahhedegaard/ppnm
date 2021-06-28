CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: uprime.pyxplot.png SIR.pyxplot.png Newton.pyxplot.png Figure8_Newton.pyxplot.png

uprime.pyxplot.png: uprime.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key bottom left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Solving \(u^{\prime\prime} = -u\)";\
	plot \
	"uprime.txt" using 1:2 with points title "\(u\)"\
	,"uprime.txt" using 1:3 with points title "\(u^\prime\)"\
	' | tee tmp.gpi | pyxplot
    

SIR.pyxplot.png: SIR.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "Time (days)";\
	set ylabel "People";\
	set title "SIR";\
	plot \
	"SIR.txt" using 1:2 with  lines title "Susceptible"\
	,"SIR.txt" using 1:3 with lines title "Infectious"\
	,"SIR.txt" using 1:4 with lines title "Recovered or dead"\
	' | tee tmp.gpi | pyxplot

Newton.pyxplot.png: Newton.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "Time";\
	set ylabel "Distance";\
	set title "Newtonian gravitational three-body problem: Path";\
	plot \
	"Newton.txt" using 1:2 with lines title  "\(x_1\)"\
	,"Newton.txt" using 1:3 with lines title "\(y_1\)"\
	,"Newton.txt" using 1:4 with lines title "\(x_2\)"\
	,"Newton.txt" using 1:5 with lines title "\(y_2\)"\
	,"Newton.txt" using 1:6 with lines title "\(x_3\)"\
	,"Newton.txt" using 1:7 with lines title "\(y_3\)"\
	' | tee tmp.gpi | pyxplot

Figure8_Newton.pyxplot.png: Newton.txt main Makefile out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "\(x_1\)";\
	set ylabel "\(y_1\)";\
	set title "Newtonian gravitational three-body problem: Reproducing fig. 8";\
	plot \
	"Newton.txt" using 2:3 with lines title "Figure 8 Solution"\
	' | tee tmp.gpi | pyxplot


out.txt: main 
	./$< > $@
    
main: functions.o main.c 
    
uprime.txt: main    

SIR.txt: main    
    
clean: 
	$(RM) main main.o functions.o *.txt *.png *.gpi
    
    
    