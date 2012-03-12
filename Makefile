
fit : main.o fit.o io.o function.o errbars.o ridge.o header.h mpfit.h
	gcc -O3 -o fit main.o fit.o io.o function.o errbars.o ridge.o -L. -lmpfit -lcfitsio -lm -llapack -lblas -ansi

main.o : main.c header.h mpfit.h
	gcc -O3 -c main.c header.h -ansi

fit.o : fit.c header.h mpfit.h
	gcc -O2 -c fit.c header.h -ansi

io.o : io.c header.h
	gcc -O2 -c io.c header.h -ansi

function.o : function.c header.h
	gcc -O3 -c function.c header.h -ansi

errbars.o : errbars.c header.h
	gcc -O3 -c errbars.c header.h -ansi

ridge.o : ridge.c header.h
	gcc -O3 -c ridge.c header.h -ansi

mpfit : mpfit.c mpfit.h
	gcc -O3 -c mpfit.c -o mpfit.o
	ar rcs libmpfit.a mpfit.o

clean : 
	rm fit *.o
