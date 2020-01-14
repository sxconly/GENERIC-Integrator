//
// Created by Xiaocheng Shang on 14/01/2020.
// Copyright © 2020 Xiaocheng Shang. All rights reserved.
//

//
// X. Shang and H. C. Öttinger: Structure-preserving integrators for dissipative systems based on reversible-irreversible splitting,
// Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, (2020).
//

//
// Damped harmonic oscillator ( U(q) = 0.5*k*q^2 ), assuming m = 1.
// Note that in order to run either of the YBABY, mYBABY, ADG, or RK3 methods, you could simply 'comment out' the other three.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

///////////////////////////////////////////////////////////////////////////

int N = 1; // number of particles

# define d 1 // dimension

# define dof d*N // degrees of freedom

# define Time 200.0 // integration time

double k_spring = 1.0; // spring constant

///////////////////////////////////////////////////////////////////////////

void Update_Force_Laplacian( double q[], double my_Force[], double my_Laplacian[] );

int main(int argc, char * argv[]){
    
/////////////
// Declare Variables
/////////////

    double d_dt = 1.30; // ratio of the increment of the stepsize
    double i_dt = 0.1/pow(d_dt, 9.0); // initial stepsize
    int number_of_stepsizes = 16;
    int number_of_runs = 1;
    
    double KbT = 1.0; // thermal energy
    double gamma = 0.01; // strength of dissipative force.
    
    int i, n;
    
    double my_Force[1], my_Laplacian[1];
            
    double Fc[dof], q[dof], p[dof], q0[dof], p0[dof], S;
    double q_ref[dof], p_ref[dof], S_ref, q0_ref[dof], p0_ref[dof];
    double q_1[dof], p_1[dof], S_1, q_2[dof], p_2[dof], S_2, q_3[dof], p_3[dof], S_3;
    
    double KE = 0.0, PE = 0.0, E0 = 0.0;
    
    double RMSE_q = 0.0, RMSE_p = 0.0, RMSE_S = 0.0, RMSE_E = 0.0;
    
    double RMSE_average_q[number_of_stepsizes];
    double RMSE_average_p[number_of_stepsizes];
    double RMSE_average_S[number_of_stepsizes];
    double RMSE_average_E[number_of_stepsizes];
    
/////////////////
// Different Stepsizes
////////////////
    
    int n_of_s; // index correspondint to 'number of stepsize'
    double dt = i_dt; // stepsize
    for (n_of_s = 0; n_of_s < number_of_stepsizes; n_of_s ++) {
        
/////////////////
// Set useful constants
////////////////
        
        // See details of the exact solutions in
        // GENERIC Integrators: Structure Preserving Time Integration for Thermodynamic Systems by H. C. Öttinger
        // https://www.degruyter.com/view/j/jnet.2018.43.issue-2/jnet-2017-0034/jnet-2017-0034.xml
        double omega = sqrt( 1.0 - 0.25*gamma*gamma ), A = cos(omega*dt), B = sin(omega*dt)/(omega*dt), D = - exp(-gamma*dt) * B*B;
        double C = ( 1.0 - exp(-gamma*dt) * ( 1.0 - 0.5*gamma*omega*sin(2.0*omega*dt) - 0.25*gamma*gamma*cos(2.0*omega*dt)) / (omega*omega) ) / (2.0*gamma*dt);
        double E = ( 1.0 - exp(-gamma*dt) * ( 1.0 + 0.5*gamma*omega*sin(2.0*omega*dt) - 0.25*gamma*gamma*cos(2.0*omega*dt)) / (omega*omega) ) / (2.0*gamma*dt*dt);
        
        double fac_1, fac_2; // Modifying factors
    
///////////////////////////////////////////////////////////////////////////
        
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
            
            KE = 0.5*p[0]*p[0], PE = 0.5*k_spring*q[0]*q[0];
            
            E0 = KE + PE + S;
            
            q_ref[0] = q[0], p_ref[0] = p[0], S_ref = S;
            
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
            
///////////////////////////////////////////////////////////////////////////
            
//////////////
// All done - let's go!
//////////////
            
            for (i = 0; i < nsteps; i++) {
                
///////////////////////////////////////////////////////////////////////////
                
                if ( i >= nsteps_per ) {
                    
                    RMSE_q = (q[0] - q_ref[0]) * (q[0] - q_ref[0]);
                    RMSE_runs_average_q += RMSE_q;
                    
                    RMSE_p = (p[0] - p_ref[0]) * (p[0] - p_ref[0]);
                    RMSE_runs_average_p += RMSE_p;
                    
                    RMSE_S = (S - S_ref) * (S - S_ref);
                    RMSE_runs_average_S += RMSE_S;
                    
                    KE = 0.5*p[0]*p[0], PE = 0.5*k_spring*q[0]*q[0];
                    
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
                
                for (n = 0; n < dof; n++) {
                    q0_ref[n] = q_ref[n];
                    p0_ref[n] = p_ref[n];
                }
                
                fac_1 = 1.0 + dt*dt*k_spring/6.0;
                fac_2 = 1.0 - dt*dt*k_spring/12.0;
                
///////////////////////////////////////////////////////////////////////////                
//////////////////////////// Exact Solution ///////////////////////////////
///////////////////////////////////////////////////////////////////////////
                
                for (n = 0; n < dof; n++) {
                    q_ref[n] = exp(-0.5*gamma*dt) * ( A*q0_ref[n] + B*( 0.5*gamma*q0_ref[n] + p0_ref[n] )*dt );
                    p_ref[n] = exp(-0.5*gamma*dt) * ( A*p0_ref[n] - B*( 0.5*gamma*p0_ref[n] + q0_ref[n] )*dt );
                    S_ref  += gamma*( C*p0_ref[n]*p0_ref[n]*dt + D*p0_ref[n]*q0_ref[n]*dt*dt + E*q0_ref[n]*q0_ref[n]*dt*dt );
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
/////////////////////////////// ADG ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
                
//                 for (n = 0; n < dof; n++) {
// 
//                     p[n] = ( 4.0*p0[n] - 4.0*dt*k_spring*q0[n] - 2.0*dt*gamma*p0[n] - dt*dt*k_spring*p0[n] ) /( 4.0 + 2.0*dt*gamma + dt*dt*k_spring );
// 
//                     q[n] = ( 4.0*q0[n] + 4.0*dt*p0[n] + 2.0*dt*gamma*q0[n] - dt*dt*k_spring*q0[n] ) /( 4.0 + 2.0*dt*gamma + dt*dt*k_spring );
// 
//                     S   += 0.25*dt*gamma*( p0[n] + p[n] )*( p0[n] + p[n] );
//                 }
                
                
///////////////////////////////////////////////////////////////////////////                
////////////////////////////// RK3 ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
                
//                 for (n = 0; n < dof; n++) {
//                     q_1[n] = dt* p0[n];
//                     p_1[n] = dt* ( - k_spring*q0[n] - gamma*p0[n] );
//                     S_1   = dt* gamma* p0[n]*p0[n]; // p0
//                 }                
//                 
//                 for (n = 0; n < dof; n++) {
//                     q_2[n] = dt* ( p0[n] + 0.5*p_1[n] );
//                     p_2[n] = dt* ( - k_spring* ( q0[n] + 0.5*q_1[n] ) - gamma* ( p0[n] + 0.5*p_1[n] ) ); 
//                     S_2   = dt* gamma* ( p0[n] + 0.5*p_1[n] ) * ( p0[n] + 0.5*p_1[n] );
//                 }
// 
//                 for (n = 0; n < dof; n++) {
//                     q_3[n] = dt* ( p0[n] - p_1[n] + 2.0*p_2[n] );
//                     p_3[n] = dt* ( - k_spring* ( q0[n] - q_1[n] + 2.0*q_2[n] ) - gamma* ( p0[n] - p_1[n] + 2.0*p_2[n] ) );
//                     S_3   = dt* gamma* ( p0[n] - p_1[n] + 2.0*p_2[n] ) * ( p0[n] - p_1[n] + 2.0*p_2[n] ); 
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
//////////////////////////// Harmonic Potential ///////////////////////////
///////////////////////////////////////////////////////////////////////////

void Update_Force_Laplacian( double q[], double my_Force[], double my_Laplacian[] )
{

/////////////////
// U(q) = 0.5*k*q^2
/////////////////

    my_Force[0] = - k_spring*q[0]; // The modulus of the Force
    my_Laplacian[0] = d*k_spring; // The Laplacian

}
