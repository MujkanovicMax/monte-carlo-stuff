#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "functions.h"



//inputs pos,sun angles, tau_c, 

void mc_1_layer(double *pos_i, double theta_s, double phi_s, double tau_c, double dz, int Ntot, double *pos_f , double *dir_f ){

    
    time_t t;
    srand((unsigned) time(&t));
    
    
    
    double dir[3];
    double pos[3];
    
    double bext = tau_c / dz;
    
    int N_up = 0;
    int N_dn = 0;
    int N_sca=0;
    int N_sca_m=0;
    
    
    double w[]={1,0,0};
    double w_alt[]={0,1,0};
    
    
    int I = 0;
    double avg = 0;
    
    while(I < Ntot){
        
        pos[0]=pos_i[0];
        pos[1]=pos_i[1];
        pos[2]=pos_i[2];
        
        dir[0]=sin(theta_s)*cos(phi_s);
        dir[1]=sin(theta_s)*sin(phi_s);
        dir[2]=-cos(theta_s);
        
        while(1==1){    
            
            ///////////////////////////Transmission through free mean path//////////////////////////////////////
            
            double tau = taufmp();
            double ds = tau/bext;
//             printf("tau=%f  ds=%f\n",tau,ds);
            double norm_dir = sqrt(dir[0]*dir[0] + dir[1]*dir[1] +dir[2]*dir[2]);
            double n = ds/norm_dir;
            
            for(int j = 0; j<3;j++){
                pos[j] = pos[j] + dir[j]*n;
            }
            
            
            if(pos[2] <= pos_i[2]-dz){
                N_dn++;
//                 printf("N_dn++ at I=%d\n\n",I);  
                
                
                
                break;
            }
            else if(pos[2] > pos_i[2]){
                N_up++;
//                 printf("N_up++ at I=%d\n\n",I);
                break;
            }
            
            ///////////////////////////Scattering in cloud layer//////////////////////////////////////////////
            
            N_sca++;
            
            norm(dir,3);
            
            
            double mu = rayscatHG();
            double phi=randphi();
//             printf("theta=%f    phi=%f \n\n",acos(mu)*180/M_PI,phi*180/M_PI);
            
            double u[3];
            double v[3];
            cross(dir,w,u);
            if(u[0] == 0 && u[1] == 0 && u[2] == 0){
                cross(dir,w_alt,u);
            }
            norm(u,3);
            cross(dir,u,v);
            
            checkorth(dir,u);
            checkorth(dir,v);
            checkorth(v,u);
            
            double dir_temp[3];
            double nphi[3];
            
            
            
            for(int k = 0;k<3;k++){
                
                nphi[k] = cos(phi)*u[k] + sin(phi)*v[k];
                
            }
            
            
            for(int k = 0;k<3;k++){
                
                dir_temp[k] = cos(M_PI/2. - acos(mu))*nphi[k] + sin(M_PI/2. - acos(mu))*dir[k];
                
            }
            //         printvec(dir,3);
            //         printf("dir\n");
            //         printvec(dir_temp,3);
            //         printf("dir_temp\n");
            
            if(checkangle(dir,dir_temp,3,mu)==1){
                
                for(int k = 0;k<3;k++){
                    dir[k]=dir_temp[k];
                }
                
            }
            else{
                printf("breaking loop at I = %d\n",I);
                break;
                
            }
            
        }   
        ///max number of scatterings
        if(N_sca>N_sca_m){
         N_sca_m = N_sca;   
        
          
        }
        ///avg number of scatterings
        avg += N_sca;
        N_sca = 0;
        ///photon count
        I++;
    }
    avg = avg/Ntot;
    printf("Trans = %f  Refl = %f   StdDev = %f   N_sca_m = %d  N_sca_avg = %f     dE = %d\n", (double) N_dn/Ntot, (double) N_up/Ntot, 
           sqrt(((double)Ntot - N_dn)/((double)Ntot*N_dn)) , N_sca_m , avg ,N_up + N_dn -I);
    printf("N_dn = %d   N_up = %d   Ntot = %d\n", N_dn, N_up, Ntot);
}