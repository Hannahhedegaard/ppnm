CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

.PHONY:default
defalut: out.error.txt out.bessel.txt out.myfun.png

out.myfun.png: out.error.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set xlabel "x";\
	set ylabel "erf(x)";\
	set title "Error function";\
	plot \
	"$<" using 1:2 with line title "error function"\
	' | pyxplot

out.error.txt: error
	./$< > $@

out.bessel.txt: bessel
	./$< > $@
