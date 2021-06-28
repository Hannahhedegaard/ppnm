CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

out.txt: hello
	./$< > $@

hello: hello.c

#hello: hello.o
#	gcc -o hello hello.o
#
#hello.o: hello.c
#	gcc -c hello.c
