
fit : main.o fit.o noise.o header.h mpfit.h libmpfit.a
	gcc -o fit main.o fit.o noise.o -L. -lmpfit -lcfitsio -lm -lgsl -lgslcblas -ansi

main.o : main.c header.h mpfit.h libmpfit.a
	gcc -c main.c header.h -lcfitsio -lm -ansi

fit.o : fit.c header.h mpfit.h libmpfit.a
	gcc -c fit.c header.h -lcfitsio -lm -ansi

noise.o : noise.c header.h
	gcc -c noise.c header.h -lm -lgsl -lgslcblas -ansi
