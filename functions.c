#include <math.h>
#include <stdio.h>
#include <stdlib.h>



double rayscat(){
    int r = rand();
    double r1 = (double)r/RAND_MAX;
    
    double q = -8*r1+4;
    double D = 1+q*q/4;
    double u = pow((-q/2+sqrt(D)),1/3.);
    
    double mur = u-1/u;
    
    return mur;
}


double rayscatHG(){
    
    int r = rand();
    double r1 = (double)r/RAND_MAX;
    double g2=0.85;
    return (pow((2*r1/(1-g2*g2)+1/(g2+1)),-2)-g2*g2-1)/(-2*g2);   
    
}

double taufmp(){
    
    int r = rand();
    double r1 = (double)r/RAND_MAX;
    
    return -log(1-r1);
    
}

double sfcrefl(){
    
    int r = rand();
    double r1 = (double)r/RAND_MAX;
    return asin(r1);
}