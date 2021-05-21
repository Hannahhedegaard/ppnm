CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: output.txt
	@echo
	@echo '=============================='
	@echo '========    TASK A    ========'
	@echo '=============================='
	@echo 'Check output below.'
	@echo
	@cat output.txt
	@echo 

output.txt: main
	./$< > $@

main: main.o integrate.o ../util.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png ../util.o

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)
