#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "functions.h"

void mc_1_layer(double *pos_i, double theta_s, double phi_s, double tau_c, double dz);


int main(){
 
    int dim = 3;
    
    double pos[]    = {0.,0.,1000.};
    
    double theta_s  = 30*M_PI/180.;
    double phi_s    = 30*M_PI/180.;
    
    double tau_c    = 10.;
    double dz       = 1000.;
    
    
    
    mc_1_layer(pos, theta_s, phi_s, tau_c, dz);
    
    
    
    
}