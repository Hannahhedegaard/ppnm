CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

#default:

#TheWay.pyxplot.png: TheWay.txt main Makefile out.txt
#	echo '\
#	set terminal png;\
#	set output "$@";\
#	set key top left;\
#	set tics out;\
#	set xlabel "x";\
#	set ylabel "y";\
#	set title "Sine function example";\
#	plot \
#	"TheWay.txt" using 1:2 with points title "u"\
#	,"TheWay.txt" using 1:3 with points title "u prime"\
#	' | tee tmp.gpi | pyxplot

out.txt: main 
	./$< > $@
    
main: functions.c main.c lineqfunctions.c
    
clean: 
	$(RM) main main.o functions.o *.txt *.png *.gpi