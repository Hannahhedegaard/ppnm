CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: Figure1.pyxplot.png Figure2.pyxplot.png

Figure1.pyxplot.png: data1.txt input1.dat main out.txt Makefile
	echo '\
	set term png;\
	set key top left;\
	set tics in;\
	set xrange [-4:4];\
	set yrange [-1.5:1.5];\
	set title "Replication of Figure 1.2 from the book";\
	set xlabel "x";\
	set ylabel "y";\
	set output "$@";\
	plot \
	"input1.dat" using 1:2 with points title "Data points" \
	,"data1.txt" using 1:2 with lines title "Akima sub-spline" \
	,"data1.txt" using 1:3 with lines title "Cubic spline" \
	' | tee tmp.gpi | pyxplot

Figure2.pyxplot.png: data2.txt input2.dat main out.txt Makefile
	echo '\
	set term png;\
	set key top right;\
	set tics in;\
	set xrange [0.5:10.5];\
	set yrange [1.25:5.25];\
	set title "Example of akima sub-spline vs. cubic spline";\
	set xlabel "x";\
	set ylabel "y";\
	set output "$@";\
	plot \
	"input2.dat" using 1:2 with points title "Data points" \
	,"data2.txt" using 1:2 with lines title "Akima sub-spline" \
	,"data2.txt" using 1:3 with lines title "Cubic spline" \
	' | tee tmp.gpi | pyxplot

out.txt: main 
	./$< > $@
      
data1.txt: main

data2.txt: main

main: functions.c main.c cubic.c
   
clean: 
	$(RM) main *.o out.txt *.png *.gpi data1.txt data2.txt