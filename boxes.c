#include <stdio.h>
#include <math.h>
#include <stdlib.h>



void checkborder(double *pos, double *dir, double *dr ,  double *border){
    
    for(int k=0;k<3;k++){
        if(dir[k] > 0){
            border[k] = ((int)pos[k]/(int)dr[k]+1)*dr[k];
        }
        else{
            border[k] = ((int)pos[k]/(int)dr[k])*dr[k];
            
        }
        
    }
    
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
    
    int x_max, y_max, z_max;
    int dx, dy, dz;
    
    x_max, y_max, z_max = 9000;
    double dr[] = {3000,3000,3000};
    
    
    
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
    
    int dim = 3;
    
    double pos[]    = {2300,1010,7021.};
    double pos_f[dim];
    double dir_f[dim];
    double dir[dim];
    double border[dim];
    
    
    double theta_s  = 60*M_PI/180.;
    double phi_s    = 60*M_PI/180.;
    dir[0]=sin(theta_s)*cos(phi_s);
    dir[1]=sin(theta_s)*sin(phi_s);
    dir[2]=-cos(theta_s);
    
    //     next box
    
    
    checkborder(pos,dir,dr,border);
    double step = choosestep(pos,dir,border);
    for(int k=0;k<3;k++){
        pos_f[k] = pos[k] + step*dir[k];
    }
    
    printf("%f  %f  %f\n", pos_f[0], pos_f[1], pos_f[2]);
    
}