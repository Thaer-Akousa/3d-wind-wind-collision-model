CC      = gcc
LDLIBS  = -lgsl -lgslcblas -lm
CFLAGS  = -Wall -g
OBJ     = numerical.o lib.o vars.o main.o
INC     = vars.h lib.h numerical.h config.h

%.o:	%.c $(INC)
	$(CC) -c $(CFLAGS) $*.c
	
bs:	$(OBJ) $(INC)
	$(CC) $(CFLAGS) -o bs $(OBJ) $(LDLIBS)

run:
	./run_bs.sh

clean:
	rm -f bs *.o *.~
