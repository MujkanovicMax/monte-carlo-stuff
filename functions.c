#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void printvec(double *v, int len){
 
    for(int i=0; i<len;i++){
     
        printf("%f  ",v[i]);
    }
    printf("\n");
}


void norm(double *v, int len){
    double n;
    for(int i=0;i<len;i++){
     
        n+=v[i]*v[i];
    }
    n=sqrt(n);
    for(int i=0;i<len;i++){
     
        v[i]=v[i]/n;
    }
}

double randphi(){
    int r = rand();
    double r1 = (double)r/RAND_MAX;
    return 2*M_PI*r1;
}

void cross(double *u, double *v, double *w){
    
    w[0]=u[1]*v[2]-u[2]*v[1];
    w[1]=u[2]*v[0]-u[0]*v[2];
    w[2]=u[0]*v[1]-u[1]*v[0];
    
}

double checkorth(double *u, double *v){
 
    if(u[0]*v[0]+u[1]*v[1]+u[2]*v[2] == 0){
     
        return 0;
        
    }
    else{
        return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
    }
    
}


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