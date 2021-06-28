CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

gamma.pyxplot.png: data_gamma.txt Makefile tabulated_gamma.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gamma function";\
	plot \
	 "$<" using 1:2 with line title "tgamma from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-gamma"\
	,"$<" using 1:4 with line title "homemade Gamma"\
	,"tabulated_gamma.txt" using 1:2 with points title "tabulated data"\
	' | pyxplot

data_gamma.txt: main
	./$< > $@

main: main.o gamma.o
