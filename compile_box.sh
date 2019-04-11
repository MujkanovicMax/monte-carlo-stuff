rm field.txt *.o boxes
CC=gcc
CFLAGS=" -Wall -g"
$CC -c $CFLAGS functions.c -lm 
$CC -c $CFLAGS randsca.c -lm
$CC -c $CFLAGS ascii.c -lm
$CC -c $CFLAGS boxes.c -lm
$CC $CFLAGS -o boxes boxes.o ascii.o randsca.o functions.o -lm

