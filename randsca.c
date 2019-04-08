#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "functions.h"





int main(){
    
    time_t t;
    srand((unsigned) time(&t));
    
    
    double pos[3];
    double dir[3];
    double dz=1000.; //=10km
    double theta_s = 15*M_PI/180.;
    double tau_c = 10.;
    double bext = tau_c / dz;
    
    int N_up = 0;
    int N_dn = 0;
    
    pos[0]=0.;
    pos[1]=0.;
    pos[2]=dz;
    
    dir[0]=sin(theta_s);
    dir[1]=0.;
    dir[2]=cos(theta_s);
    double w[]={1,0,0};
    
    for(int i=0;i<1;i++){
        
     
        double tau = taufmp();
        double ds = tau/bext;
        double norm_dir = sqrt(dir[0]*dir[0] + dir[1]*dir[1] +dir[2]*dir[2]);
        double n = ds/norm_dir;
        for(int j = 0; j<3;j++){
            pos[i] = pos[i] + dir[i]*n;
        }
        
        
        if(pos[2] <= 0){
            N_dn++;
        }
        else if(pos[2] >= dz){
            N_up++;
        }
        
        for(int k=0;k<3;k++){
            dir[k]=dir[k]/norm_dir;
        }
        
        
        double mu = rayscatHG();
        double u[3];
        double v[3];
        cross(dir,w,u);
        cross(dir,u,v);
        double dir_temp[3];
        for(int k=0; k<3;k++){
            dir_temp[k] = dir[k] + u[k]*tan(acos(mu));
        }
        norm(dir_temp,3);
        
        double phi=randphi();
        
        for(int k=0; k<3;k++){
            dir[k] = dir_temp[k] + tan(phi)*v[k];
        }
        norm(dir,3);
        printvec(dir,3);
        printf("%f\n",dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
        
//         printf("%f  %f  %f\n",checkorth(u,v) ,checkorth(u,dir), checkorth(dir,v));
        
        
        
        
    }
    
}