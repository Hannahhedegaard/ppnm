CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


erf.pyxplot.png: data_erf.txt Makefile tabulated_erf.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "homemade Erf"\
	,"tabulated_erf.txt" using 1:2 with points title "tabulated data"\
	' | pyxplot

data_erf.txt: main
	./$< > $@

main: main.o erf.o
