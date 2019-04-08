rm *.o randsca
CC=gcc
CFLAGS=" -Wall"
$CC -c $CFLAGS functions.c -lm
$CC -c $CFLAGS randsca.c -lm
$CC $CFLAGS -o randsca randsca.o functions.o -lm

