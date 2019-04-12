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

void hor_boundary(double *pos, double *bounds){
    
    for(int k = 0;k<2;k++){
        if(pos[k]>=bounds[k]-0.0001){
            pos[k] = 0;   
        }
        else if(pos[k]<=0+0.0001){
            pos[k] = bounds[k];
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


void checkbox(double *pos, double *dir, double *dr, int *lmn){
    for(int k=0;k<3;k++){
        lmn[k] = (int)(floor(pos[k]/dr[k]));
    }
    
}

double choosestep(double *pos, double *dir, double *border,int *coo){
    
    double vstep[3];
    double tmp;
    for(int k=0;k<3;k++){
        if(dir[k]==0){
            vstep[k]=100000000000; //pro progaming
            
        }
        else{
            vstep[k] = (border[k]-pos[k])/dir[k];
        }
    }
    tmp=vstep[0];
    
    for(int k = 0;k<3;k++){
        if(tmp > vstep[k]){
            tmp = vstep[k];
            *coo = k;
        }
        
    }
    
    return tmp;
    
}

int n_betacheck(double *pos, double *z,int nlvl){
    int tmp;
    for(int k =0;k<nlvl;k++){
        if(pos[2]>=z[k]){
            tmp=k;
            break;
        }
        else if(pos[2]<1){
            
            tmp = 49;
        }
    }
    return tmp;
}

/*
 * i *n*t cloud(int *lmn){
 * if(lmn[0]>=20 && lmn[0]<=60 && lmn[1]>=20 && lmn[1]<=60 && lmn[2]>=5 && lmn[2]<=20){
 *    return 1;
 *    }
 *    else{
 *        return 0;
 *        }
 *        }
 */


int main(){
    
    //     set box par
    
    time_t t;
    srand((unsigned) time(&t));
    
    double *z;
    double *k_s;
    double *k_a;
    int nlvl;
    
    read_3c_file ("./ksca_kabs/ksca_kabs_lambda350.dat", &z, &k_s, &k_a, &nlvl);
    
    double x_step=1;
    double y_step=1;
    double z_step=1;
    
    double x = nlvl*x_step;
    double y = nlvl*y_step;
    
    int nlyr = nlvl-1;
    //     double beta_a[nlyr];
    double beta_s[nlyr];
    double beta_ext[nlyr];
    double tau_lyr;
    double w0[nlyr];
    
    
    double dx=x_step;
    double dy=y_step;
    double dz=z_step;
    
    
    
    
    
    
    
    
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
    
    int N_sca_up=0;
    int N_abs=0;
    int N_dn[120][120];
    int N_up=0;
    int Ntot=100000;
    
    double pos[dim];
    double bounds[] = {120,120,120};
    double pos_f[dim];
    //     double dir_f[dim];
    double dir[dim];
    double border[dim];
    double beta_ext_lyr;
    double beta_s_lyr;
    double beta_s_c=20;
    double w0_lyr;
    
    
    for(int k=0;k<120;k++){
        for(int h=0;h<120;h++){
            N_dn[k][h] = 0;
        }
    }
    
    
    
    //     printf("%f  %f  %f\n",pos_f[0],pos_f[1],pos_f[2]);
    
    
    //     photon loop
    int I=0;
    int N_sca_m=0;
    int coo=0;
    int n_beta;
    
    
    
    
    
    while(I<Ntot){
        
        
        
        
        double tau = taufmp();
        int lmn[3];
        double path [3];
        
        int N_sca=0;
        
        pos[0]=randnum()*120;
        pos[1]=randnum()*120;
        pos[2]=120;
        
        double theta_s  = 15*M_PI/180.;
        double phi_s    = 0*M_PI/180.;
        
        dir[0]=sin(theta_s)*cos(phi_s);
        dir[1]=sin(theta_s)*sin(phi_s);
        dir[2]=-cos(theta_s);
        
        double dr[dim];
        int T = 0;
        
        
        
        checkbox(pos,dir,dr,lmn);
        
        
        dr[0]=dx;
        dr[1]=dy;
        dr[2]=dz;
        
        
        while(1==1){
            n_beta=n_betacheck(pos,z,nlvl);
            beta_s_lyr=k_s[n_beta];
            beta_ext_lyr= beta_s_lyr+k_a[n_beta];
            
            
            //             if(cloud(lmn)==1){
            T=0;
            if(pos[1]>30){
                if(pos[1]<80){
                    if(pos[0]>30){
                        if(pos[0]<80){
                            if(pos[2]>30){
                                if(pos[2]<80){
                                    T=1;
                                }
                            }
                        }
                    }
                }
            }
            
            if(T==1){        
                
                beta_ext_lyr=beta_ext_lyr + beta_s_c;
                beta_s_lyr=beta_s_lyr + beta_s_c;
            }
            
            w0_lyr=beta_s_lyr/beta_ext_lyr;
            
            checkborder(pos,dir,dr,border);
            double step = choosestep(pos,dir,border,&coo);
            //             printf("step= %f\n",step);
            for(int k=0;k<3;k++){
                path[k]=dir[k]*step;
                pos_f[k] = pos[k] + path[k];
            }
            //             printvec(pos_f,3);
            double ds = givelen(path,3);
            for(int q = 0;q<3;q++){
                
                //                 printf("%d\n",lmn[q]);
                
            }
            //             printvec(pos,3);
            //             printvec(pos_f,3);
            
            tau_lyr = beta_ext_lyr*ds;
            if(tau_lyr>= tau){
                double tmp_vec[3];
                double b = randnum();
                
                if(b<w0_lyr){                    
                    step = tau/tau_lyr*step;
                    for(int k=0;k<3;k++){
                        path[k]=dir[k]*step;
                        pos[k] = pos[k] + path[k];
                    }
                    
                    N_sca++;
                    //                     if(cloud(lmn)==1){
                    if(T==1){
                        double c=randnum();
                        if(c<beta_s_c/beta_s_lyr){
                            scattering(pos,dir,tmp_vec);
                            //                             printf("I=%d  ",I);
                            //                             printvec(tmp_vec,3);
                        }
                        else{
                            scattering_ray(pos,dir,tmp_vec);
                            
                        }
                        tau=taufmp();
                    }
                    else{
                        scattering_ray(pos,dir,tmp_vec);
                        tau = taufmp();
                        
                        for(int k=0;k<3;k++){
                            dir[k] = tmp_vec[k];
                        }
                    }
                }
                else{
                    N_abs++;
                    break;
                }
                
                
            }
            
            else{
                
                tau=tau-tau_lyr;
                
                if(pos_f[2]>=z[0]-0.1){
                    N_up++;
                    
                    break;
                }
                else if(pos_f[2]<=0.0001){
                    double a =randnum();
                    double A=0.13;
                    if(a<A){
                        
                        double theta_sfc = sfcrefl();
                        double phi_sfc = randphi();
                        
                        dir[0]=sin(theta_sfc)*cos(phi_sfc);
                        dir[1]=sin(theta_sfc)*sin(phi_sfc);
                        dir[2]=cos(theta_sfc);
                    }
                    else{
                        checkbox(pos,dir,dr,lmn);
                        N_dn[lmn[0]][lmn[1]]++;
                        break;
                    }
                }
                
                
                
                for(int k =0;k<3;k++){
                    pos[k]=pos_f[k];
                }
                //                 printvec(pos,3);
                hor_boundary(pos,bounds);
                checkbox(pos,dir,dr,lmn);
                
                
            }            
            if(N_sca>N_sca_m){
                N_sca_m = N_sca;   
            }
            
        }
        I++;
        //         if(I%1000==0){
        //             printf("I=%d",I);
        //         }
    }
    //     printf("N_dn = %d   N_abs = %d N_up = %d  I = %d  N_sca_up = %d\n", N_dn, N_abs, N_up, I, N_sca_up);
    FILE *f;
    f = fopen("field.txt","w+");
    
    for(int k=0;k<120;k++){
        fprintf(f,"\n");
        for(int h=0;h<120;h++){
            fprintf(f,"%d ", N_dn[k][h]); 
        }
    }
    fclose(f);
}