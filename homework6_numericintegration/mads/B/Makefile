CC = gcc-10
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: output.txt
	@echo
	@echo '=============================='
	@echo '========    TASK B    ========'
	@echo '=============================='
	@echo 'Check output below. Note the dramatic decrease in the number of function calls when using'
	@echo 'the Clenshaw-Curtis transformation for integrable singularities.'
	@echo 'For the integral of 4*sqrt(1 - x*x) we can basically get pi up to machine precision (double).'
	@echo 'For this integral the CC transformation uses more evaluation but gives higher precision.'
	@echo
	@cat output.txt
	@echo 

output.txt: main
	./$< > $@

main: main.o ccintegrate.o ../A/integrate.o ../util.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png ../util.o

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
