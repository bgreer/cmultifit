
fit : main.o fit.o noise.o io.o header.h mpfit.h
	gcc -o fit main.o fit.o noise.o io.o -L. -lmpfit -lcfitsio -lm -lgsl -lgslcblas -ansi

main.o : main.c header.h mpfit.h
	gcc -c main.c header.h -ansi

fit.o : fit.c header.h mpfit.h
	gcc -c fit.c header.h -ansi

noise.o : noise.c header.h
	gcc -c noise.c header.h -ansi

io.o : io.c header.h
	gcc -c io.c header.h -ansi

