#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "ascii.h"
#include "functions.h"

void checkborder(double *pos, double *dir, double *dr ,  double *border){
    
    for(int k=0;k<3;k++){
        if(dir[k] > 0){
            border[k] = (floor(pos[k]/dr[k])+1)*dr[k];
        }
        else{
            if (floor(pos[k]/dr[k]) == pos[k]/dr[k] ){
                border[k] = (floor(pos[k]/dr[k])-1)*dr[k];
            }
            else{
                border[k] = (floor(pos[k]/dr[k]))*dr[k];
            }
        }
        
    }
    
}

double checklyr(double *pos, double *dir, double dz){
    
    double z_border=0;
    
    if(dir[2] > 0){
        z_border = (floor(pos[2]/dz)+1)*dz;
    }
    else{
        if (floor(pos[2]/dz) == pos[2]/dz ){
            z_border = (floor(pos[2]/dz)-1)*dz;
        }
        else{
            z_border = (floor(pos[2]/dz))*dz;
        }
    }
    return z_border;
}

double choosestep(double *pos, double *dir, double *border){
    
    double vstep[3];
    double tmp;
    for(int k=0;k<3;k++){
        
        vstep[k] = (border[k]-pos[k])/dir[k];
        
    }
    
    for(int k = 0;k<2;k++){
        
        if(vstep[k]<vstep[k+1]){
            tmp = vstep[k];
        }
        else{
            
            tmp = vstep[k+1];
        }
        
    }
    
    return tmp;
    
}



int main(){
    
    //     set box par
    
    time_t t;
    srand((unsigned) time(&t));
    
    double *z;
    double *k_s;
    double *k_a;
    int nlvl;
    
    read_3c_file ("./ksca_kabs/ksca_kabs_lambda300.dat", &z, &k_s, &k_a, &nlvl);
    
    
    int nlyr = nlvl-1;
    double beta_a[nlyr];
    double beta_s[nlyr];
    double beta_ext[nlyr];
    double tau_lyr;
    double w0[nlyr];
    
    double dz[nlyr];
    
    
    for(int i = 0; i<nlyr;i++){
        
        dz[i] = z[i]-z[i+1];
        beta_ext[i] = k_a[i+1] + k_s[i+1];
        
    }
    
    for(int i = 0; i<nlyr;i++){
        beta_a[i] = k_a[i+1];
        beta_s[i] = k_s[i+1];
        w0[i]     = beta_s[i]/beta_ext[i];
    }
    
    
    
    
    
    
    
    
    
    
    //     box def properties
    //     double tau_def = 1.;
    //     
    //     double tau_box[(int)(x_max/dx)][(int)(y_max/dy)][(int)(z_max/dz)];
    //     for(int i = 0;i<x_max/dx;i++){
    //         for(int j = 0;j<y_max/dy;j++){
    //             for(int k = 0;k<z_max/dz;k++){
    //                 
    //                 tau_box[i][j][k] = tau_def;
    //             }
    //         }
    //     }
    //     //      box special props
    //     
    //     tau_box[1][1][1] = 100.;
    //     
    //      initialization params
    
    /////////////////////////////////////////////////
    
    int dim = 3;
    
    int N_abs=0;
    int N_dn=0;
    int N_up=0;
    int Ntot=10;
    
    double pos[]    = {1,10,120};
    double pos_f[dim];
    double dir_f[dim];
    double dir[dim];
    double border;
    
    
    
    
    
    //     printf("%f  %f  %f\n",pos_f[0],pos_f[1],pos_f[2]);
    
    
    //     photon loop
    int I=0;
    int N_sca_m=0;
    
    
    while(I<Ntot){
        
        
        
        double tau = taufmp();
        int j =0;
        double path [3];
        
        int N_sca=0;
        
        pos[0]=1;
        pos[1]=10;
        pos[2]=120;
        
        double theta_s  = 15*M_PI/180.;
        double phi_s    = 15*M_PI/180.;
        
        dir[0]=sin(theta_s)*cos(phi_s);
        dir[1]=sin(theta_s)*sin(phi_s);
        dir[2]=-cos(theta_s);
        
        
        while(1==1){
            
//             printvec(pos,3);
            
            double z_border = checklyr(pos,dir,dz[j]);
            double step = fabs(dz[j]/dir[2]);
            //         printf("step = %f\n", step);
            
            for(int k=0;k<3;k++){
                path[k]=dir[k]*step;
                pos_f[k] = pos[k] + path[k];
            }
            //         printvec(path,3);
            double ds = givelen(path,3);
            //         printf("ds=%f\n",ds);
            
            tau_lyr = beta_ext[j]*ds;
            //         printf("tau_lyr=%f tau=%f\n",tau_lyr,tau);
            if(tau_lyr > tau){
                //             printf("hello");
                double tmp_vec[3];
                double a = randnum();
                
                if(a<w0[j]){
                    printf("a = %f  w0 = %f j=%d\n",a,w0[j],j);
                    N_sca++;
                    scattering_ray(pos,dir,tmp_vec);
                    //                 printf("scattering");
                    
                    
                    for(int k=0;k<3;k++){
                        dir[k] = tmp_vec[k];
                    }            
                    step = tau/tau_lyr*step;
                    for(int k=0;k<3;k++){
                        path[k]=dir[k]*step;
                        pos[k] = pos[k] + path[k];
                    }
                    
                }
                else{
                    N_abs++;
                    break;
                }
                
                
            }
            
            else{
                tau=tau-tau_lyr;
                if(pos_f[2]>=z[0]-0.0001){
                    N_up++;
                    
                    break;
                }
                else if(pos_f[2]<=0.0001){
                    
                    N_dn++;
                    
                    break;
                    
                }
                
                for(int k =0;k<3;k++){
                    pos[k]=pos_f[k];
                }
                
                //         lyr update
                //             printvec(pos_f,3);
                
                if (dir[2]<0){
                    j++;
                }
                else{
                    j--;
                }
                //             printf("j=%d\n",j);
            }
            
            if(N_sca>N_sca_m){
                N_sca_m = N_sca;   
            }
            
        }
        I++;
        
    }
    printf("N_dn = %d   N_up = %d  I = %d  N_sca_m = %d\n", N_dn, N_up,I, N_sca_m);
}