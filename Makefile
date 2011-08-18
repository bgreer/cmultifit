
fit : main.o fit.o noise.o io.o function.o header.h mpfit.h
	gcc -O3 -o fit main.o fit.o noise.o io.o function.o -L. -lmpfit -lcfitsio -lm -lgsl -lgslcblas -ansi

main.o : main.c header.h mpfit.h
	gcc -O3 -c main.c header.h -ansi

fit.o : fit.c header.h mpfit.h
	gcc -O3 -c fit.c header.h -ansi

noise.o : noise.c header.h
	gcc -O3 -c noise.c header.h -ansi

io.o : io.c header.h
	gcc -O3 -c io.c header.h -ansi

function.o : function.c header.h
	gcc -c function.c header.h -ansi
