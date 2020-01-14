%%
%% Created by Xiaocheng Shang on 14/01/2020.
%% Copyright © 2020 Xiaocheng Shang. All rights reserved.
%%

%%
%% X. Shang and H. C. Öttinger: Structure-preserving integrators for dissipative systems based on reversible-irreversible splitting,
%% Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, (2020).
%%

%%
%% Damped Nonlinear Oscillator ( U(q) = - k*cos(c*q) )
%% solved by using the average discrete gradient (ADG) method
%%

clc;
clear all;
close all;
 
%% 
d_dt = 1.3; % ratio of the increment of the stepsize
number_of_stepsizes = 16;
RMSE = zeros(number_of_stepsizes,5);
 
dt = 0.1/d_dt^9;

%% Parameters 
k = 1; c = 1; m = 1; T = 1; gamma = 0.01;
     
for j=1:number_of_stepsizes
     
    Time = 180;
    nsteps = floor(Time/dt); % number of steps
    sols = zeros(nsteps,7);
    
    %% Initail conditions
    sols(1,1) = 2; % q
    sols(1,2) = 0; % p
    sols(1,3) = 0; % S
    sols(1,4) = 0.5*sols(1,2)^2/m; % KE
    sols(1,5) = - k*cos( c*sols(1,1) ); % PE
    sols(1,6) = sols(1,4) + sols(1,5); % KE + PE
    sols(1,7) = sols(1,6) + T*sols(1,3); % E
     
    E0 = sols(1,7); % initial total energy
     
    %% Reference solutions
     
    refs = zeros(nsteps,7);
    refs(1,:) = sols(1,:);
     
    %%
    for i=1:nsteps
        
        %% Initial values at the beginning of each step
        q0 = sols(i,1);
        p0 = sols(i,2);
		
	    %% The Newton-Raphson method
		
        % Set up the iteration
        tolerence = 1.e8;
        xx = [0.5, 0.5]; % initial guesses
        iter = 0;
        maxIterations = 40;
        
        % Begin the iteration
        while tolerence > 1.e-12
            
            iter = iter + 1;
            
            q = xx(1);
            p = xx(2);
            
            % Calculate the functions
            f(1) = q - q0 - 0.5*dt/m * ( p + p0 );
            f(2) = p - p0 - dt*k * ( cos(c*q) - cos(c*q0) ) / (q - q0) + 0.5*dt*gamma * ( p + p0 );
            
            % Calculate the Jacobian
            J(1,1) = 1;
            J(1,2) = - 0.5*dt/m;
            J(2,1) = dt*k* ( c*sin(c*q)*(q - q0) + cos(c*q) - cos(c*q0) ) / (q - q0)^2;
            J(2,2) = 1 + 0.5*dt*gamma;
            
            % Solve the linear equations
            u = -J\f';
            
            % Update the solution
            xx = xx + u';
            
            % Calculate the norm
            tolerence = sqrt( u(1)*u(1) + u(2)*u(2) );

            if (iter > maxIterations)
                s = sprintf('****Did not converge within %3.0f iterations.****', maxIterations);
                disp(s)
            end
        end
        
        %% Update the solution	
        sols(i+1,1) = xx(1); % q
        sols(i+1,2) = xx(2); % p
        sols(i+1,3) = sols(i,3) + dt*gamma/(4*m*T) * ( sols(i,2) + sols(i+1,2) )^2; % S
        sols(i+1,4) = 0.5*sols(i+1,2)^2/m; % KE
        sols(i+1,5) = - k*cos( c*sols(i+1,1) ); % PE
        sols(i+1,6) = sols(i+1,4) + sols(i+1,5); % KE + PE
        sols(i+1,7) = sols(i+1,6) + T*sols(i+1,3); % E
         
    end
     
    %% Reference solutions
     
    load GENERIC_Int_Nonlinear_d001_gamma_d01_k1_N1_T180_q0-2_RK3.mat % load -> x
     
    h_ref = 0.001;
     
    for i=1:nsteps
         
        N_index = round(i*dt/h_ref);
         
        refs(i+1,1) = x(N_index,1); % q
        refs(i+1,2) = x(N_index,2); % p
        refs(i+1,3) = x(N_index,3); % S
        refs(i+1,4) = x(N_index,4); % KE
        refs(i+1,5) = x(N_index,5); % PE
        refs(i+1,6) = x(N_index,6); % KE + PE
        refs(i+1,7) = x(N_index,7); % E
    end            
     
    %% RMSE
     
    RMSE(j,1) = dt;
    RMSE(j,2) = sqrt( mse( sols(:,1) - refs(:,1) ) ); % q
    RMSE(j,3) = sqrt( mse( sols(:,2) - refs(:,2) ) ); % p
    RMSE(j,4) = sqrt( mse( sols(:,3) - refs(:,3) ) ); % S
    RMSE(j,5) = sqrt( mse( sols(:,7) - E0 ) ); % E
     
    dt = dt*d_dt;
     
end
