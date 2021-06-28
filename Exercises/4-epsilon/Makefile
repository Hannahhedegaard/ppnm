CC = clang
CFLAGS = -Wall -fwrapv
LDFLAGS = -I.
LDLIBS = -lm

out.txt: epsilon equal
	./epsilon > out.txt
	./equal >> out.txt

epsilon.o: epsilon.c 
	$(CC) $(CFLAGS) -c epsilon.c -o epsilon.o

epsilon: epsilon.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

equal.o: equal.c
	$(CC) $(CFLAGS) -c equal.c -o equal.o

equal: equal.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

