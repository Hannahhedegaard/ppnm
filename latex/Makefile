CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

ms.pdf: ms.tex fig-pyxplot.pdf
	pdflatex $<
	pdflatex $<
#	bibtex $<
#	pdflatex $<

out.data.txt:main
	./$< > $@

main:main.o mygamma.o

fig-pyxplot.pdf: fig-pyxplot.gpi out.data.txt
	pyxplot $<
