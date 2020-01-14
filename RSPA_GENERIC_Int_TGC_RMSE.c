//
// Created by Xiaocheng Shang on 12/01/2020.
// Copyright © 2020 Xiaocheng Shang. All rights reserved.
//

//
// X. Shang and H. C. Öttinger: Structure-preserving integrators for dissipative systems based on reversible-irreversible splitting,
// Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, (2020).
//

//
// Two gas containers exchanging heat and volume ( E = 0.5*p^2/m + E1 + E2 ), assuming m = 1, N*k_B = 1.
// Note that in order to run either of the YBABY, mYBABY, or RK3 methods, you could simply 'comment out' the other two.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

///////////////////////////////////////////////////////////////////////////

int N = 1; // number of particles

# define density 0.84 // particle density

# define d 1 // dimension

# define dof d*N // degrees of freedom

# define Time 30.0 // integration time

///////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]){
    
/////////////
// Declare Variables
/////////////
    
    double d_dt = 1.30; // ratio of the increment of the stepsize
    double i_dt = 0.1/pow(d_dt, 9.0); // initial stepsize
    int number_of_stepsizes = 10; // to compare the effect of using different stepsizes
    int number_of_runs = 1;
    
    int i, n;
            
    double q[dof], p[dof], q0[dof], p0[dof], S1[dof], S2[dof], E1[dof], E2[dof], E1_m[dof], E2_m[dof], A1[dof], A2[dof];
    double q_ref[dof], p_ref[dof], S_ref, q0_ref[dof], p0_ref[dof];

    double q_1[dof], p_1[dof], q_2[dof], p_2[dof], q_3[dof], p_3[dof];
    double S1_0[dof], S2_0[dof], S1_1[dof], S2_1[dof], S1_2[dof], S2_2[dof], S1_3[dof], S2_3[dof];
    double E1_0[dof], E2_0[dof], A1_0[dof], A2_0[dof], E1_1[dof], E2_1[dof], A1_1[dof], A2_1[dof], E1_2[dof], E2_2[dof], E1_3[dof], E2_3[dof];
    
    double KE = 0.0, PE = 0.0, E0 = 0.0;
    
    double RMSE_q = 0.0, RMSE_p = 0.0, RMSE_S = 0.0, RMSE_E = 0.0;
    
    double RMSE_average_q[number_of_stepsizes];
    double RMSE_average_p[number_of_stepsizes];
    double RMSE_average_S1[number_of_stepsizes];
    double RMSE_average_S2[number_of_stepsizes];
    double RMSE_average_S[number_of_stepsizes];
    double RMSE_average_E[number_of_stepsizes];
    
///////////////////////////////////////////////////////////////////////////
    
        double h_ref = 0.001;
        int N_steps = (int) floor(Time/h_ref) + 1;
        double q_i[N_steps], p_i[N_steps], S_i[N_steps], E_i[N_steps];
        
///////////////////////////////////////////////////////////////////////////
        // Import initial configurations from data files
        FILE *fp_q;
        fp_q = fopen("GENERIC_Int_TGC_d001_N1_T30_L1_p2_E4_ad5_RK3_q.txt", "r");
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
        fp_p = fopen("GENERIC_Int_TGC_d001_N1_T30_L1_p2_E4_ad5_RK3_p.txt", "r");
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
        fp_S = fopen("GENERIC_Int_TGC_d001_N1_T30_L1_p2_E4_ad5_RK3_S.txt", "r");
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

        fp_E = fopen("GENERIC_Int_TGC_d001_N1_T30_L1_p2_E4_ad5_RK3_E.txt", "r");
        if(fp_E == NULL){
            printf("Open file failure!");
            exit(1);
        }
        for (n = 0; n < N_steps; n++) {
            fscanf(fp_E, "%lf,", &E_i[n]); // IMPORTANT, "%lf" (if the output is of type double), NOT "%f"!
        }
        fclose(fp_E); // NOT "close"!!!
        
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
        RMSE_average_S1[n_of_s] = 0.0;
        RMSE_average_S2[n_of_s] = 0.0;
        RMSE_average_S[n_of_s] = 0.0;
        RMSE_average_E[n_of_s] = 0.0;
        
        // number of steps
        int nsteps = (int) floor(Time/dt) + 1;
        
        double per = 1.0; int nsteps_per = (int) floor((1.0-per)*nsteps); // only collect last 100*per % data
        
        double RMSE_runs_average_q = 0.0;
        double RMSE_runs_average_p = 0.0;
        double RMSE_runs_average_S1 = 0.0;
        double RMSE_runs_average_S2 = 0.0;
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
            
            // N=1/1.381*10^(23); h=6.626*10^(-34); c=exp(1)^(2.5)/N*(4*pi/(3*N*h^2))^(1.5)=2.544*10^(44) = e^102.2476703501216; ==> Const = 102.2476703501216;
            // N=1; h=6.626*10^(-34); c=exp(1)^(2.5)/N*(4*pi/(3*N*h^2))^(1.5) = 3.590*10^(101) = e^233.8392935112113; ==> Const = 233.8392935112113;
            double Ac = 1.0, Lg = 1.0, Const = 102.2476703501216, alpha = 0.5;
            A1[0] = 0.0, A2[0] = 0.0; // Modifying factors: alpha_1 and alpha_2;
            q[0] = 1.0, p[0] = 2.0, E1[0] = 2.0, E2[0] = 2.0;
            
            KE = 0.5*p[0]*p[0];
            E0 = KE + E1[0] + E2[0];
            
            S1[0] = 1.5 * log( E1[0] ) + log( q[0] * Ac ) + Const;
            S2[0] = 1.5 * log( E2[0] ) + log( ( 2.0*Lg - q[0] ) * Ac) + Const;
            
            double q_ref[nsteps], p_ref[nsteps], S_ref[nsteps], E_ref[nsteps];
    
            for (n = 0; n < nsteps; n++) {
                
                int N_index = (int) round(n*dt/h_ref);
                
                q_ref[n] = q_i[N_index];
                p_ref[n] = p_i[N_index];
                S_ref[n] = S_i[N_index];
                E_ref[n] = E_i[N_index];
            }
            
//////////////
// All done - let's go!
//////////////

            for (i = 0; i < nsteps; i++) {
                
///////////////////////////////////////////////////////////////////////////
                
                if ( i >= nsteps_per ) {
                    
                    // Important to notice the index of the reference solution
                    RMSE_q = (q[0] - q_ref[i]) * (q[0] - q_ref[i]);
                    RMSE_runs_average_q += RMSE_q;
                    
                    RMSE_p = (p[0] - p_ref[i]) * (p[0] - p_ref[i]);
                    RMSE_runs_average_p += RMSE_p;
                    
                    RMSE_S = (S1[0] + S2[0] - S_ref[i]) * (S1[0] + S2[0] - S_ref[i]);
                    RMSE_runs_average_S += RMSE_S;

                    
                    E1_m[0] = exp( (2.0/3.0) * ( S1[0] - log( q[0] * Ac ) - Const ) );
                    E2_m[0] = exp( (2.0/3.0) * ( S2[0] - log( ( 2.0*Lg - q[0] ) * Ac) - Const ) );
                    
                    KE = 0.5*p[0]*p[0];

                    RMSE_E = (KE + E1_m[0] + E2_m[0] - E0) * (KE + E1_m[0] + E2_m[0] - E0);
                    RMSE_runs_average_E += RMSE_E;
                    
                }
                
///////////////////////////////////////////////////////////////////////////
////////////////////////////// MAIN LOOP //////////////////////////////////
///////////////////////////////////////////////////////////////////////////

                for (n = 0; n < dof; n++) {
                    q0[n] = q[n];
                    p0[n] = p[n];
                    S1_0[n] = S1[n];
                    S2_0[n] = S2[n];
                }
                
///////////////////////////////////////////////////////////////////////////              
/////////////////////////////// YBABY /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
                
                // Midpoint method for half a step
                for (n = 0; n < dof; n++) {
                    q0[n] = q[n];
                    p0[n] = p[n];
                    S1_0[n] = S1[n];
                    S2_0[n] = S2[n];
                }
                
                for (n = 0; n < dof; n++) {
                    E1_0[n] = exp( (2.0/3.0) * ( S1_0[n] - log( q0[n] * Ac ) - Const ) );
                    E2_0[n] = exp( (2.0/3.0) * ( S2_0[n] - log( ( 2.0*Lg - q0[n] ) * Ac) - Const ) );
                    
                    S1_1[n] = S1_0[n] + 0.25*dt* 9.0 * alpha / (4.0*E1_0[n]) * ( 1.0/E1_0[n] - 1.0/E2_0[n] );
                    S2_1[n] = S2_0[n] - 0.25*dt* 9.0 * alpha / (4.0*E2_0[n]) * ( 1.0/E1_0[n] - 1.0/E2_0[n] );
                }
                
                for (n = 0; n < dof; n++) {
                    
                    E1_1[n] = exp( (2.0/3.0) * ( S1_1[n] - log( q0[n] * Ac ) - Const ) );
                    E2_1[n] = exp( (2.0/3.0) * ( S2_1[n] - log( ( 2.0*Lg - q0[n] ) * Ac) - Const ) );
                    
                    S1[n] = S1_0[n] + 0.5*dt* 9.0 * alpha / (4.0*E1_1[n]) * ( 1.0/E1_1[n] - 1.0/E2_1[n] );
                    S2[n] = S2_0[n] - 0.5*dt* 9.0 * alpha / (4.0*E2_1[n]) * ( 1.0/E1_1[n] - 1.0/E2_1[n] );
                }
                
                // Verlet method for a step
                for (n = 0; n < dof; n++) {
                    
                    E1[n] = exp( (2.0/3.0) * ( S1[n] - log( q[n] * Ac ) - Const ) );
                    E2[n] = exp( (2.0/3.0) * ( S2[n] - log( ( 2.0*Lg - q[n] ) * Ac) - Const ) );
                    
                    // p
                    p[n] += dt* ( E1[n] / q[n] - E2[n] / ( 2.0*Lg - q[n] ) ) / 3.0;
                    
                    // q
                    q[n] += dt* p[n];
                    
                    E1[n] = exp( (2.0/3.0) * ( S1[n] - log( q[n] * Ac ) - Const ) );
                    E2[n] = exp( (2.0/3.0) * ( S2[n] - log( ( 2.0*Lg - q[n] ) * Ac) - Const ) );
                    
                    // p
                    p[n] += dt* ( E1[n] / q[n] - E2[n] / ( 2.0*Lg - q[n] ) ) / 3.0;
                    
                }
                
                // Midpoint method for half a step
                for (n = 0; n < dof; n++) {
                    q0[n] = q[n];
                    p0[n] = p[n];
                    S1_0[n] = S1[n];
                    S2_0[n] = S2[n];
                }
                
                for (n = 0; n < dof; n++) {
                    E1_0[n] = exp( (2.0/3.0) * ( S1_0[n] - log( q0[n] * Ac ) - Const ) );
                    E2_0[n] = exp( (2.0/3.0) * ( S2_0[n] - log( ( 2.0*Lg - q0[n] ) * Ac) - Const ) );
                    
                    S1_1[n] = S1_0[n] + 0.25*dt* 9.0 * alpha / (4.0*E1_0[n]) * ( 1.0/E1_0[n] - 1.0/E2_0[n] );
                    S2_1[n] = S2_0[n] - 0.25*dt* 9.0 * alpha / (4.0*E2_0[n]) * ( 1.0/E1_0[n] - 1.0/E2_0[n] );
                }
                
                for (n = 0; n < dof; n++) {
                    
                    E1_1[n] = exp( (2.0/3.0) * ( S1_1[n] - log( q0[n] * Ac ) - Const ) ); // q_1 -> q0
                    E2_1[n] = exp( (2.0/3.0) * ( S2_1[n] - log( ( 2.0*Lg - q0[n] ) * Ac) - Const ) ); // q_1 -> q0
                    
                    S1[n] = S1_0[n] + 0.5*dt* 9.0 * alpha / (4.0*E1_1[n]) * ( 1.0/E1_1[n] - 1.0/E2_1[n] );
                    S2[n] = S2_0[n] - 0.5*dt* 9.0 * alpha / (4.0*E2_1[n]) * ( 1.0/E1_1[n] - 1.0/E2_1[n] );
                }
                       
   
///////////////////////////////////////////////////////////////////////////                
/////////////////////////////// mYBABY ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
                
//                 // Midpoint method for half a step
//                 for (n = 0; n < dof; n++) {
//                     q0[n] = q[n];
//                     p0[n] = p[n];
//                     S1_0[n] = S1[n];
//                     S2_0[n] = S2[n];
//                 }
//                     
//                 for (n = 0; n < dof; n++) {
//                     E1_0[n] = exp( (2.0/3.0) * ( S1_0[n] - log( q0[n] * Ac ) - Const ) );
//                     E2_0[n] = exp( (2.0/3.0) * ( S2_0[n] - log( ( 2.0*Lg - q0[n] ) * Ac) - Const ) );
//                     
//                     A1_0[n] = 1.0 + dt*dt* ( 5.0*p0[n]*p0[n] / q0[n]              - 2.0*( E1_0[n] / q0[n] - E2_0[n] / ( 2.0*Lg - q0[n] ) ) ) / ( 54.0*q0[n] );
//                     A2_0[n] = 1.0 + dt*dt* ( 5.0*p0[n]*p0[n] / ( 2.0*Lg - q0[n] ) + 2.0*( E1_0[n] / q0[n] - E2_0[n] / ( 2.0*Lg - q0[n] ) ) ) / ( 54.0* ( 2.0*Lg - q0[n] ) );
//                     
//                     S1_1[n] = S1_0[n] + 0.25*dt* 9.0 * alpha * A2_0[n] / (4.0*E1_0[n]) * ( A2_0[n]/E1_0[n] - A1_0[n]/E2_0[n] );
//                     S2_1[n] = S2_0[n] - 0.25*dt* 9.0 * alpha * A1_0[n] / (4.0*E2_0[n]) * ( A2_0[n]/E1_0[n] - A1_0[n]/E2_0[n] );
//                 }
//                 
//                 for (n = 0; n < dof; n++) {
//                     
//                     E1_1[n] = exp( (2.0/3.0) * ( S1_1[n] - log( q0[n] * Ac ) - Const ) ); // q_1 -> q0
//                     E2_1[n] = exp( (2.0/3.0) * ( S2_1[n] - log( ( 2.0*Lg - q0[n] ) * Ac) - Const ) ); // q_1 -> q0
//                     
//                     A1_1[n] = 1.0 + dt*dt* ( 5.0*p0[n]*p0[n] / q0[n]              - 2.0*( E1_1[n] / q0[n] - E2_1[n] / ( 2.0*Lg - q0[n] ) ) ) / ( 54.0*q0[n] );
//                     A2_1[n] = 1.0 + dt*dt* ( 5.0*p0[n]*p0[n] / ( 2.0*Lg - q0[n] ) + 2.0*( E1_1[n] / q0[n] - E2_1[n] / ( 2.0*Lg - q0[n] ) ) ) / ( 54.0* ( 2.0*Lg - q0[n] ) );
//                     
//                     S1[n] = S1_0[n] + 0.5*dt* 9.0 * alpha * A2_1[n] / (4.0*E1_1[n]) * ( A2_1[n]/E1_1[n] - A1_1[n]/E2_1[n] );
//                     S2[n] = S2_0[n] - 0.5*dt* 9.0 * alpha * A1_1[n] / (4.0*E2_1[n]) * ( A2_1[n]/E1_1[n] - A1_1[n]/E2_1[n] );
//                 }
//                 
//                 // Verlet method for a step
//                 for (n = 0; n < dof; n++) {
//                     
//                     E1[n] = exp( (2.0/3.0) * ( S1[n] - log( q[n] * Ac ) - Const ) );
//                     E2[n] = exp( (2.0/3.0) * ( S2[n] - log( ( 2.0*Lg - q[n] ) * Ac) - Const ) );
//                     
//                     // p
//                     p[n] += dt* ( E1[n] / q[n] - E2[n] / ( 2.0*Lg - q[n] ) ) / 3.0;
//                     
//                     // q
//                     q[n] += dt* p[n];
//                     
//                     E1[n] = exp( (2.0/3.0) * ( S1[n] - log( q[n] * Ac ) - Const ) );
//                     E2[n] = exp( (2.0/3.0) * ( S2[n] - log( ( 2.0*Lg - q[n] ) * Ac) - Const ) );
//                     
//                     // p
//                     p[n] += dt* ( E1[n] / q[n] - E2[n] / ( 2.0*Lg - q[n] ) ) / 3.0;
//                     
//                 }
//                 
//                 // Midpoint method for half a step
//                 for (n = 0; n < dof; n++) {
//                     q0[n] = q[n];
//                     p0[n] = p[n];
//                     S1_0[n] = S1[n];
//                     S2_0[n] = S2[n];
//                 }
//                     
//                 for (n = 0; n < dof; n++) {
//                     E1_0[n] = exp( (2.0/3.0) * ( S1_0[n] - log( q0[n] * Ac ) - Const ) );
//                     E2_0[n] = exp( (2.0/3.0) * ( S2_0[n] - log( ( 2.0*Lg - q0[n] ) * Ac) - Const ) );
//                     
//                     A1_0[n] = 1.0 + dt*dt* ( 5.0*p0[n]*p0[n] / q0[n]              - 2.0*( E1_0[n] / q0[n] - E2_0[n] / ( 2.0*Lg - q0[n] ) ) ) / ( 54.0*q0[n] );
//                     A2_0[n] = 1.0 + dt*dt* ( 5.0*p0[n]*p0[n] / ( 2.0*Lg - q0[n] ) + 2.0*( E1_0[n] / q0[n] - E2_0[n] / ( 2.0*Lg - q0[n] ) ) ) / ( 54.0* ( 2.0*Lg - q0[n] ) );
//                     
//                     S1_1[n] = S1_0[n] + 0.25*dt* 9.0 * alpha * A2_0[n] / (4.0*E1_0[n]) * ( A2_0[n]/E1_0[n] - A1_0[n]/E2_0[n] );
//                     S2_1[n] = S2_0[n] - 0.25*dt* 9.0 * alpha * A1_0[n] / (4.0*E2_0[n]) * ( A2_0[n]/E1_0[n] - A1_0[n]/E2_0[n] );
//                 }
//                 
//                 for (n = 0; n < dof; n++) {
//                     
//                     E1_1[n] = exp( (2.0/3.0) * ( S1_1[n] - log( q0[n] * Ac ) - Const ) ); // q_1 -> q0
//                     E2_1[n] = exp( (2.0/3.0) * ( S2_1[n] - log( ( 2.0*Lg - q0[n] ) * Ac) - Const ) ); // q_1 -> q0
//                     
//                     A1_1[n] = 1.0 + dt*dt* ( 5.0*p0[n]*p0[n] / q0[n]              - 2.0*( E1_1[n] / q0[n] - E2_1[n] / ( 2.0*Lg - q0[n] ) ) ) / ( 54.0*q0[n] );
//                     A2_1[n] = 1.0 + dt*dt* ( 5.0*p0[n]*p0[n] / ( 2.0*Lg - q0[n] ) + 2.0*( E1_1[n] / q0[n] - E2_1[n] / ( 2.0*Lg - q0[n] ) ) ) / ( 54.0* ( 2.0*Lg - q0[n] ) );
//                     
//                     S1[n] = S1_0[n] + 0.5*dt* 9.0 * alpha * A2_1[n] / (4.0*E1_1[n]) * ( A2_1[n]/E1_1[n] - A1_1[n]/E2_1[n] );
//                     S2[n] = S2_0[n] - 0.5*dt* 9.0 * alpha * A1_1[n] / (4.0*E2_1[n]) * ( A2_1[n]/E1_1[n] - A1_1[n]/E2_1[n] );
//                 }
                
 
///////////////////////////////////////////////////////////////////////////                
////////////////////////////// RK3 ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
                
//                 for (n = 0; n < dof; n++) {
// 
//                     E1_0[n] = exp( (2.0/3.0) * ( S1_0[n] - log( q0[n] * Ac ) - Const ) );
//                     E2_0[n] = exp( (2.0/3.0) * ( S2_0[n] - log( ( 2.0*Lg - q0[n] ) * Ac) - Const ) );
// 
//                     q_1[n]  =   dt* p0[n];
//                     p_1[n]  =   dt* (2.0/3.0) * ( E1_0[n] / q0[n] - E2_0[n] / ( 2.0*Lg - q0[n] ) );
//                     S1_1[n] =   dt* 9.0 * alpha / (4.0*E1_0[n]) * ( 1.0/E1_0[n] - 1.0/E2_0[n] );
//                     S2_1[n] = - dt* 9.0 * alpha / (4.0*E2_0[n]) * ( 1.0/E1_0[n] - 1.0/E2_0[n] );
//                 }
// 
//                 for (n = 0; n < dof; n++) {
// 
//                     E1_1[n] = exp( (2.0/3.0) * ( ( S1_0[n] + 0.5*S1_1[n] ) - log( ( q0[n] + 0.5*q_1[n] ) * Ac ) - Const ) );
//                     E2_1[n] = exp( (2.0/3.0) * ( ( S2_0[n] + 0.5*S2_1[n] ) - log( ( 2.0*Lg - ( q0[n] + 0.5*q_1[n] ) ) * Ac) - Const ) );
// 
//                     q_2[n]  =   dt* ( p0[n] + 0.5*p_1[n] );
//                     p_2[n]  =   dt* (2.0/3.0) * ( E1_1[n] / ( q0[n] + 0.5*q_1[n] ) - E2_1[n] / ( 2.0*Lg - ( q0[n] + 0.5*q_1[n] ) ) );
//                     S1_2[n] =   dt* 9.0 * alpha / (4.0*E1_1[n]) * ( 1.0/E1_1[n] - 1.0/E2_1[n] );
//                     S2_2[n] = - dt* 9.0 * alpha / (4.0*E2_1[n]) * ( 1.0/E1_1[n] - 1.0/E2_1[n] );
//                 }
// 
//                 for (n = 0; n < dof; n++) {
// 
//                     E1_2[n] = exp( (2.0/3.0) * ( ( S1_0[n] - S1_1[n] + 2.0*S1_2[n] ) - log( ( q0[n] - q_1[n] + 2.0*q_2[n] ) * Ac ) - Const ) );
//                     E2_2[n] = exp( (2.0/3.0) * ( ( S2_0[n] - S2_1[n] + 2.0*S2_2[n] ) - log( ( 2.0*Lg - ( q0[n] - q_1[n] + 2.0*q_2[n] ) ) * Ac) - Const ) );
// 
//                     q_3[n]  =   dt* ( p0[n] - p_1[n] + 2.0*p_2[n] );
//                     p_3[n]  =   dt* (2.0/3.0) * ( E1_2[n] / ( q0[n] - q_1[n] + 2.0*q_2[n] ) - E2_2[n] / ( 2.0*Lg - ( q0[n] - q_1[n] + 2.0*q_2[n] ) ) );
//                     S1_3[n] =   dt* 9.0 * alpha / (4.0*E1_2[n]) * ( 1.0/E1_2[n] - 1.0/E2_2[n] );
//                     S2_3[n] = - dt* 9.0 * alpha / (4.0*E2_2[n]) * ( 1.0/E1_2[n] - 1.0/E2_2[n] );
//                 }
// 
//                 for (n = 0; n < dof; n++) {
//                     q[n]  += ( q_1[n] + 4.0*q_2[n] + q_3[n] ) / 6.0;
//                     p[n]  += ( p_1[n] + 4.0*p_2[n] + p_3[n] ) / 6.0;
//                     S1[n] += ( S1_1[n] + 4.0*S1_2[n] + S1_3[n] ) / 6.0;
//                     S2[n] += ( S2_1[n] + 4.0*S2_2[n] + S2_3[n] ) / 6.0;
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

        // Print the averaged quantities!
        printf("%f, %f, %f, %f, %f \n", dt, RMSE_average_q[n_of_s], RMSE_average_p[n_of_s], RMSE_average_S[n_of_s], RMSE_average_E[n_of_s]);
        
///////////////////////////////////////////////////////////////////////////
        
        dt *= d_dt;
        
    }
    
// Finished!
    return 0;
    
}
