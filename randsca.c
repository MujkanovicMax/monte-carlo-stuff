#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

int main(){
    
    time_t t;
    
    srand((unsigned) time(&t));
    
    int r = rand();
    double r1 = r/RAND_MAX;
//     double phir = 2*M_PI *r;
    
    printf("%f    %d\n",r1, r);
    
}