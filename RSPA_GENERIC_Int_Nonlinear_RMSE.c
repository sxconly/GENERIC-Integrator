//
// Created by Xiaocheng Shang on 12/01/2020.
// Copyright © 2020 Xiaocheng Shang. All rights reserved.
//

//
// X. Shang and H. C. Öttinger: Structure-preserving integrators for dissipative systems based on reversible-irreversible splitting,
// Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, (2020).
//

//
// Damped nonlinear oscillator ( U(q) = - k*cos(c*q) ), assuming m = 1.
// Note that in order to run either of the YBABY, mYBABY, or RK3 methods, you could simply 'comment out' the other two.
//
        
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

///////////////////////////////////////////////////////////////////////////

int N = 1; // number of particles

# define d 1 // dimension

# define dof d*N // degrees of freedom

# define Time 180.0 // integration time

double k_spring = 1.0; // spring constant

///////////////////////////////////////////////////////////////////////////

void Update_Force_Laplacian( double q[], double my_Force[], double my_Laplacian[] );

int main(int argc, char * argv[]){
    
/////////////
// Declare Variables
/////////////
    
    double d_dt = 1.30; // ratio of the increment of the stepsize
    double i_dt = 0.1/pow(d_dt, 9.0); // initial stepsize
    int number_of_stepsizes = 16; // to compare the effect of using different stepsizes
    int number_of_runs = 1;
    
    double KbT = 1.0; // thermal energy
    double gamma = 0.01; // strength of dissipative force.
    
    int i, n;
    
    double my_Force[1], my_Laplacian[1];
            
    double Fc[dof], q[dof], p[dof], q0[dof], p0[dof], S, fac_1, fac_2;
    double q_ref[dof], p_ref[dof], S_ref, q0_ref[dof], p0_ref[dof];
    double q_1[dof], p_1[dof], S_1, q_2[dof], p_2[dof], S_2, q_3[dof], p_3[dof], S_3;
            
    double KE = 0.0, PE = 0.0, E0 = 0.0;
    
    double RMSE_q = 0.0, RMSE_p = 0.0, RMSE_S = 0.0, RMSE_E = 0.0;
    
    double RMSE_average_q[number_of_stepsizes];
    double RMSE_average_p[number_of_stepsizes];
    double RMSE_average_S[number_of_stepsizes];
    double RMSE_average_E[number_of_stepsizes];
    
///////////////////////////////////////////////////////////////////////////
    
        double h_ref = 0.001;
        int N_steps = (int) floor(Time/h_ref) + 1;
        double q_i[N_steps], p_i[N_steps], S_i[N_steps], E_i[N_steps];
        
///////////////////////////////////////////////////////////////////////////
        // Import initial configurations from data files
        FILE *fp_q;
        fp_q = fopen("GENERIC_Int_Nonlinear_d001_gamma_d01_k1_N1_T180_q0-2_RK3_q.txt", "r");
        if(fp_q == NULL){
            printf("Open file failure!");
            exit(1);
        }
        for (n = 0; n < N_steps; n++) {
            fscanf(fp_q, "%lf,", &q_i[n]); // IMPORTANT, "%lf" (if the output is of type double), NOT "%f"!
        }
        fclose(fp_q);
        
///////////////////////////////////////////////////////////////////////////
        // Import initial momentum from data files
        FILE *fp_p;
        fp_p = fopen("GENERIC_Int_Nonlinear_d001_gamma_d01_k1_N1_T180_q0-2_RK3_p.txt", "r");
        if(fp_p == NULL){
            printf("Open file failure!");
            exit(1);
        }
        for (n = 0; n < N_steps; n++) {
            fscanf(fp_p, "%lf,", &p_i[n]); // IMPORTANT, "%lf" (if the output is of type double), NOT "%f"!
        }
        fclose(fp_p);
        
///////////////////////////////////////////////////////////////////////////
        // Import initial entropy from data files
        FILE *fp_S;
        fp_S = fopen("GENERIC_Int_Nonlinear_d001_gamma_d01_k1_N1_T180_q0-2_RK3_S.txt", "r");
        if(fp_S == NULL){
            printf("Open file failure!");
            exit(1);
        }
        for (n = 0; n < N_steps; n++) {
            fscanf(fp_S, "%lf,", &S_i[n]); // IMPORTANT, "%lf" (if the output is of type double), NOT "%f"!
        }
        fclose(fp_S);
        
///////////////////////////////////////////////////////////////////////////
        // Import initial entropy from data files
        FILE *fp_E;
        fp_E = fopen("GENERIC_Int_Nonlinear_d001_gamma_d01_k1_N1_T180_q0-2_RK3_E.txt", "r");
        if(fp_E == NULL){
            printf("Open file failure!");
            exit(1);
        }
        for (n = 0; n < N_steps; n++) {
            fscanf(fp_E, "%lf,", &E_i[n]); // IMPORTANT, "%lf" (if the output is of type double), NOT "%f"!
        }
        fclose(fp_E);
        
///////////////////////////////////////////////////////////////////////////

    
/////////////////
// Different Stepsizes
////////////////
    
    int n_of_s; // index correspondint to 'number of stepsize'
    double dt = i_dt; // stepsize
    for (n_of_s = 0; n_of_s < number_of_stepsizes; n_of_s ++) {
        
/////////////////
// Set useful constants
////////////////
        
        // very important!
        RMSE_average_q[n_of_s] = 0.0;
        RMSE_average_p[n_of_s] = 0.0;
        RMSE_average_S[n_of_s] = 0.0;
        RMSE_average_E[n_of_s] = 0.0;
        
        // very important!        
        RMSE_average_q[n_of_s] = 0.0;
        RMSE_average_p[n_of_s] = 0.0;
        RMSE_average_S[n_of_s] = 0.0;
        RMSE_average_E[n_of_s] = 0.0;
        
        
        // number of steps
        int nsteps = (int) floor(Time/dt) + 1;
        
        double per = 1.0; int nsteps_per = (int) floor((1.0-per)*nsteps); // only collect last 100*per % data
        
        double RMSE_runs_average_q = 0.0;
        double RMSE_runs_average_p = 0.0;
        double RMSE_runs_average_S = 0.0;
        double RMSE_runs_average_E = 0.0;
        
/////////////////
// Average a couple of runs (if needed)
////////////////
        
        int n_of_r; // index correspondint to 'number of run'
        for (n_of_r = 0; n_of_r < number_of_runs; n_of_r++) {

////////////////
// Set Initial Conditions
////////////////
            
            q[0] = 2.0, p[0] = 0.0, S = 0.0;
            KE = 0.5*p[0]*p[0], PE = - k_spring*cos(q[0]);
            E0 = KE + PE + S;

            double q_ref[nsteps], p_ref[nsteps], S_ref[nsteps], E_ref[nsteps];
            
            for (n = 0; n < nsteps; n++) {
                
                int N_index = (int) round(n*dt/h_ref);
                
                q_ref[n] = q_i[N_index];
                p_ref[n] = p_i[N_index];
                S_ref[n] = S_i[N_index];
                E_ref[n] = E_i[N_index];
            }
            
////////////////
// Calculate Initial Forces
////////////////
            
            for (n = 0; n < dof; n++) {
                Fc[n] = 0.0;
            }
            
            Update_Force_Laplacian( q, my_Force, my_Laplacian );
            
            for (n = 0; n < dof; n++) {
                Fc[n] += my_Force[n];
            }
            
//////////////
// All done - let's go!
//////////////
            
            for (i = 0; i < nsteps; i++) {
                
                if ( i >= nsteps_per ) {
                    
                    // Important to notice the index of the reference solution
                    RMSE_q = (q[0] - q_ref[i]) * (q[0] - q_ref[i]);
                    RMSE_runs_average_q += RMSE_q;
                    
                    RMSE_p = (p[0] - p_ref[i]) * (p[0] - p_ref[i]);
                    RMSE_runs_average_p += RMSE_p;
                    
                    RMSE_S = (S - S_ref[i]) * (S - S_ref[i]);
                    RMSE_runs_average_S += RMSE_S;
                    
                    KE = 0.5*p[0]*p[0], PE = - k_spring*cos(q[0]);

                    RMSE_E = (KE + PE + KbT*S - E0) * (KE + PE + KbT*S - E0);
                    RMSE_runs_average_E += RMSE_E;
                    
                }
   
///////////////////////////////////////////////////////////////////////////
////////////////////////////// MAIN LOOP //////////////////////////////////
///////////////////////////////////////////////////////////////////////////                
                
                for (n = 0; n < dof; n++) {
                    q0[n] = q[n];
                    p0[n] = p[n];
                }
                
///////////////////////////////////////////////////////////////////////////             
/////////////////////////////// YBABY /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
                
                for (n = 0; n < dof; n++) {

                    // Exact solution for half a step
                    p0[n] = p[n]; // Important!
                    p[n]  = exp( - 0.5*gamma*dt )*p0[n];
                    S    += 0.5*p0[n]*p0[n] * ( 1.0 - exp(-gamma*dt ) ); // Note p0[n]!
                    
                    // Verlet method for a step
                    p[n] += 0.5*dt* my_Force[0];
                    
                    q[n] += dt* p[n]; 
                    Update_Force_Laplacian( q, my_Force, my_Laplacian ); 

                    p[n] += 0.5*dt* my_Force[0];

                    // Exact solution for half a step
                    p0[n] = p[n]; // Important!
                    p[n]  = exp( - 0.5*gamma*dt )*p0[n];
                    S    += 0.5*p0[n]*p0[n] * ( 1.0 - exp(-gamma*dt ) ); // Note p0[n]!
                    
                }
                
                
///////////////////////////////////////////////////////////////////////////                
/////////////////////////////// mYBABY ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
                
//                 for (n = 0; n < dof; n++) {
// 
//                     // Exact solution for half a step
//                     p0[n] = p[n]; // Important!
//                     p[n]  = exp( - 0.5*gamma*dt * fac_1 )*p0[n];
//                     S    += 0.5*p0[n]*p0[n] * fac_1 * ( 1.0 - exp(-gamma*dt*fac_1 ) ); // Note p0[n]!
//                     
//                     // Verlet method for a step
//                     p[n] += 0.5*dt* my_Force[0];
//                     
//                     q[n] += dt* p[n]; 
//                     Update_Force_Laplacian( q, my_Force, my_Laplacian ); 
//                     fac_1 = 1.0 + dt*dt*my_Laplacian[0]/6.0;
// 
//                     p[n] += 0.5*dt* my_Force[0];
// 
//                     // Exact solution for half a step
//                     p0[n] = p[n]; // Important!
//                     p[n]  = exp( - 0.5*gamma*dt * fac_1 )*p0[n];
//                     S    += 0.5*p0[n]*p0[n] * fac_1 * ( 1.0 - exp(-gamma*dt*fac_1 ) ); // Note p0[n]!
//                     
//                 }
            
                
///////////////////////////////////////////////////////////////////////////                
////////////////////////////// RK3 ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
                
//                 for (n = 0; n < dof; n++) {
//                     q_1[n] = dt* p0[n];
//                     p_1[n] = dt* ( - k_spring*sin(q0[n]) - gamma*p0[n] );
//                     S_1    = dt* gamma* p0[n]*p0[n]; // p0
//                 }
// 
//                 for (n = 0; n < dof; n++) {
//                     q_2[n] = dt* ( p0[n] + 0.5*p_1[n] );
//                     p_2[n] = dt* ( - k_spring* sin( q0[n] + 0.5*q_1[n] ) - gamma* ( p0[n] + 0.5*p_1[n] ) );
//                     S_2    = dt* gamma* ( p0[n] + 0.5*p_1[n] ) * ( p0[n] + 0.5*p_1[n] );
//                 }
// 
//                 for (n = 0; n < dof; n++) {
//                     q_3[n] = dt* ( p0[n] - p_1[n] + 2.0*p_2[n] );
//                     p_3[n] = dt* ( - k_spring* sin( q0[n] - q_1[n] + 2.0*q_2[n] ) - gamma* ( p0[n] - p_1[n] + 2.0*p_2[n] ) );
//                     S_3    = dt* gamma* ( p0[n] - p_1[n] + 2.0*p_2[n] ) * ( p0[n] - p_1[n] + 2.0*p_2[n] );
//                 }
// 
//                 for (n = 0; n < dof; n++) {
//                     q[n] += ( q_1[n] + 4.0*q_2[n] + q_3[n] ) / 6.0;
//                     p[n] += ( p_1[n] + 4.0*p_2[n] + p_3[n] ) / 6.0;
//                     S    += ( S_1 + 4.0*S_2 + S_3 ) / 6.0;
//                 }

                
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
                
            }

        }
        
        RMSE_runs_average_q /= number_of_runs;
        RMSE_runs_average_p /= number_of_runs;
        RMSE_runs_average_S /= number_of_runs;
        RMSE_runs_average_E /= number_of_runs;
        
///////////////////////////////////////////////////////////////////////////

        RMSE_average_q[n_of_s] = sqrt( RMSE_runs_average_q / (nsteps - nsteps_per) );
        RMSE_average_p[n_of_s] = sqrt( RMSE_runs_average_p / (nsteps - nsteps_per) );
        RMSE_average_S[n_of_s] = sqrt( RMSE_runs_average_S / (nsteps - nsteps_per) );
        RMSE_average_E[n_of_s] = sqrt( RMSE_runs_average_E / (nsteps - nsteps_per) );
        
///////////////////////////////////////////////////////////////////////////
        
        // Print the averaged quantities!
        printf("%f, %f, %f, %f, %f \n", dt, RMSE_average_q[n_of_s], RMSE_average_p[n_of_s], RMSE_average_S[n_of_s], RMSE_average_E[n_of_s]);

///////////////////////////////////////////////////////////////////////////
        
        dt *= d_dt;
        
    }
    
// Finished!
    return 0;
    
}

///////////////////////////////////////////////////////////////////////////
//////////////////////////// Nonlinear Potential //////////////////////////
///////////////////////////////////////////////////////////////////////////

void Update_Force_Laplacian( double q[], double my_Force[], double my_Laplacian[] )
{
    
/////////////////
// U(q) = - k*cos(q)
/////////////////
    
    my_Force[0] = - k_spring*sin(q[0]); // The modulus of the Force!
    my_Laplacian[0] = k_spring*cos(q[0]); // The Laplacian
    
}
