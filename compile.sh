rm *.o mcsim
CC=gcc
CFLAGS=" -Wall"
$CC -c $CFLAGS functions.c -lm
$CC -c $CFLAGS randsca.c -lm
$CC -c $CFLAGS mcsim.c -lm
$CC $CFLAGS -o mcsim mcsim.o randsca.o functions.o -lm

