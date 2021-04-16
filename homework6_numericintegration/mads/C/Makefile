CC = gcc-10
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: output.txt
	@echo
	@echo '=============================='
	@echo '========    TASK C    ========'
	@echo '=============================='
	@echo 'Check output below.'
	@echo
	@cat output.txt
	@echo 

output.txt: main
	./$< > $@

main: main.o infintegrate.o ../util.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png ../util.o

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
