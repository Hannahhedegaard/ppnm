CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

out.txt: main 
	./$< > $@
    
main: functions.c main.c 
    
clean: 
	$(RM) main *.txt *.png *.gpi