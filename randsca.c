#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

int main(){
    
    time_t t;
    int len=100000;
    
    double g=0.;
    double g2=0.85;
    double mur[len], murHG[len], taur[len], thetar[len];
    srand((unsigned) time(&t));
    
    for(int i=0;i<len;i++){
        int r = rand();
        double r1 = (double)r/RAND_MAX;
        
        double q = -8*r+4;
        double D = 1+q*q/4;
        double u = pow((-q/2+sqrt(D)),1/3.);
        
        mur[i] = u-1/u;
        murHG[i] = (pow((2*r/(1-g2*g2)+1/(g2+1)),-2)-g2*g2-1)/(-2*g2);
        taur[i] = -log(1-r);
        thetar[i] = asin(r);
    }
    
}