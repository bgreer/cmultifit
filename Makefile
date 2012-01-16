
fit : main.o fit.o io.o function.o header.h mpfit.h
	gcc -O2 -o fit main.o fit.o io.o function.o -L. -lmpfit -lcfitsio -lm -lgsl -lgslcblas -ansi

main.o : main.c header.h mpfit.h
	gcc -O2 -c main.c header.h -ansi

fit.o : fit.c header.h mpfit.h
	gcc -O2 -c fit.c header.h -ansi

io.o : io.c header.h
	gcc -O2 -c io.c header.h -ansi

function.o : function.c header.h
	gcc -c function.c header.h -ansi
