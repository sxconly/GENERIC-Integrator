%%
%% Created by Xiaocheng Shang on 14/01/2020.
%% Copyright © 2020 Xiaocheng Shang. All rights reserved.
%%

%%
%% X. Shang and H. C. Öttinger: Structure-preserving integrators for dissipative systems based on reversible-irreversible splitting,
%% Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, (2020).
%%

%%
%% Two gas containers exchanging heat and volume ( E = 0.5*p^2/m + E1 + E2 )
%% solved by using the average discrete gradient (ADG) method
%%

clc;
clear all;
close all;

%%
d_dt = 1.3; % ratio of the increment of the stepsize
number_of_stepsizes = 10;
RMSE = zeros(number_of_stepsizes,5);
 
dt = 0.1/d_dt^9;

%% Parameters
%% N = 1/1.381*10^(23); h = 6.626*10^(-34); 
%% C = exp(1)^(2.5)/N*(4*pi/(3*N*h^2))^(1.5) = 2.544*10^(44) = exp(102.2476703501216);
Ac = 1; Lg = 1; alpha = 0.5; m = 1; NkB = 1; C = exp(102.2476703501216); B = (C*Ac)^(-2/3); 
 
for j=1:number_of_stepsizes
     
    Time = 30;
    nsteps = floor(Time/dt); % number of steps
    sols = zeros(nsteps,7);
    
    %% Initail conditions
    sols(1,1) = 1; % q
    sols(1,2) = 2; % p
 
    sols(1,6) = 2; % E1
    sols(1,7) = 2; % E2
    sols(1,8) = 0.5*sols(1,2)^2/m + sols(1,6) + sols(1,7); % KE + E1 + E2
     
    E0 = sols(1,8); % initial total energy
     
    sols(1,3) = NkB*log( C *          sols(1,1)   * Ac * sols(1,6)^(1.5) ); % S1
    sols(1,4) = NkB*log( C * ( 2*Lg - sols(1,1) ) * Ac * sols(1,7)^(1.5) ); % S2
    sols(1,5) = sols(1,3) + sols(1,4); % S1 + S2
     
    %% Reference solutions
     
    refs = zeros(nsteps,8);
    refs(1,:) = sols(1,:);
     
    %%
    for i=1:nsteps
	
	    %% Initial values at the beginning of each step
	    q_0  = sols(i,1);
        p_0  = sols(i,2);
        S1_0 = sols(i,3);
        S2_0 = sols(i,4);        
		
	    %% The Newton-Raphson method
		
        % Set up the iteration
        tolerence = 1.e8;
        xx = [1, 2, 105, 105]; % initial guesses
        iter = 0;
        maxIterations = 100;
		
        % Begin the iteration
        while tolerence > 1.e-12
		
            iter = iter + 1;
			
            q  = xx(1);
            p  = xx(2);
            S1 = xx(3);
            S2 = xx(4);
			
            % Calculate the functions
            f(1) = q - q_0 - 0.5 * dt / m * ( p + p_0 );
            f(2) = p - p_0 + dt * B / 3 * (  ( 2*Lg - q_0 )^(-5/3) * exp(2*S2_0/3/NkB) + (2*Lg - q)^(-5/3) * exp(2*S2/3/NkB) - q_0^(-5/3) * exp(2*S1_0/3/NkB) - q^(-5/3) * exp(2*S1/3/NkB)  );
            f(3) = S1 - S1_0 - dt * alpha * ( 3*NkB/(2*B) )^2 * (  ( (q_0+q)/2 )^(4/3) * exp(-2*(S1_0+S1)/3/NkB) - ( Lg*(q_0+q) - (q_0+q)^2/4 )^(2/3) * exp(-(S1_0+S1+S2_0+S2)/3/NkB)  );
            f(4) = S2 - S2_0 + dt * alpha * ( 3*NkB/(2*B) )^2 * (  ( Lg*(q_0+q) - (q_0+q)^2/4 )^(2/3) * exp(-(S1_0+S1+S2_0+S2)/3/NkB) - ( 2*Lg - (q_0+q)/2 )^(4/3) * exp(-2*(S2_0+S2)/3/NkB)  );
            
			% Calculate the Jacobian
            J(1,1) = 1;
            J(1,2) = - 0.5 * dt / m;
            J(1,3) = 0;
            J(1,4) = 0;
            J(2,1) = 5 * dt * B / 9 * (  (2*Lg - q)^(-8/3) * exp(2*S2/3/NkB) + q^(-8/3) * exp(2*S1/3/NkB)  );
            J(2,2) = 1;
            J(2,3) = - dt  *B / 3 * q^(-5/3) * (2/3/NkB) * exp(2*S1/3/NkB);
            J(2,4) =   dt * B / 3 * (2*Lg - q)^(-5/3) * (2/3/NkB) * exp(2*S2/3/NkB);
            J(3,1) = - dt * alpha * ( 3*NkB/(2*B) )^2 * (2/3) * (  ( (q_0+q)/2 )^(1/3) * exp(-2*(S1_0+S1)/3/NkB) - ( Lg*(q_0+q) - (q_0+q)^2/4 )^(-1/3) * ( Lg - (q_0+q)/2 ) * exp(-(S1_0+S1+S2_0+S2)/3/NkB)  );
            J(3,2) = 0;
            J(3,3) = 1 + dt * alpha * ( 3*NkB/(4*B^2) ) * (  2*( (q_0+q)/2 )^(4/3) * exp(-2*(S1_0+S1)/3/NkB) - ( Lg*(q_0+q) - (q_0+q)^2/4 )^(2/3) * exp(-(S1_0+S1+S2_0+S2)/3/NkB)  );
            J(3,4) =   - dt * alpha * ( 3*NkB/(4*B^2) ) * ( Lg*(q_0+q) - (q_0+q)^2/4 )^(2/3) * exp(-(S1_0+S1+S2_0+S2)/3/NkB);
            J(4,1) = dt * alpha * ( 3*NkB/(2*B) )^2 * (2/3) * (  ( Lg*(q_0+q) - (q_0+q)^2/4 )^(-1/3) * ( Lg - (q_0+q)/2 ) * exp(-(S1_0+S1+S2_0+S2)/3/NkB) + ( 2*Lg - (q_0+q)/2 )^(1/3) * exp(-2*(S2_0+S2)/3/NkB)  );
            J(4,2) = 0;
            J(4,3) =   - dt * alpha * ( 3*NkB/(4*B^2) ) * ( Lg*(q_0+q) - (q_0+q)^2/4 )^(2/3) * exp(-(S1_0+S1+S2_0+S2)/3/NkB);
            J(4,4) = 1 + dt * alpha * ( 3*NkB/(4*B^2) ) * (  2*( 2*Lg - (q_0+q)/2 )^(4/3) * exp(-2*(S2_0+S2)/3/NkB) - ( Lg*(q_0+q) - (q_0+q)^2/4 )^(2/3) * exp(-(S1_0+S1+S2_0+S2)/3/NkB)  );
            
			% Solve the linear equations
            u = -J\f'; 
			
            % Update the solution
            xx = xx + u';
			
            % Calculate the norm
            tolerence = sqrt( u(1)*u(1) + u(2)*u(2) + u(3)*u(3) + u(4)*u(4) );

            if (iter > maxIterations)
                s = sprintf('****Did not converge within %3.0f iterations.****', maxIterations);
                disp(s)
            end
        end
        
        %% Update the solution		
        sols(i+1,1) = xx(1); % q
        sols(i+1,2) = xx(2); % p
        sols(i+1,3) = xx(3); % S1
        sols(i+1,4) = xx(4); % S2
        sols(i+1,5) = sols(i+1,3) + sols(i+1,4); % S1 + S2
        sols(i+1,6) = ( exp( sols(i+1,3)/NkB ) / ( C *          sols(i+1,1)   * Ac ) )^(2/3); % E1
        sols(i+1,7) = ( exp( sols(i+1,4)/NkB ) / ( C * ( 2*Lg - sols(i+1,1) ) * Ac ) )^(2/3); % E2
        sols(i+1,8) = 0.5*sols(i+1,2)^2/m + sols(i+1,6) + sols(i+1,7); % KT + E1 + E2
         
    end
     
    %% Reference solutions
    
    load GENERIC_Int_TGC_d001_N1_T30_L1_p2_E4_ad5_RK3.mat % load -> x
    
    h_ref = 0.001;
    
    for i=1:nsteps
        
        N_index = round(i*dt/h_ref);
        
        refs(i+1,1) = x(N_index,1); % q
        refs(i+1,2) = x(N_index,2); % p
        refs(i+1,3) = x(N_index,3); % S1
        refs(i+1,4) = x(N_index,4); % S2
        refs(i+1,5) = x(N_index,5); % S1 + S2
        refs(i+1,6) = x(N_index,6); % E1
        refs(i+1,7) = x(N_index,7); % E2
        refs(i+1,8) = x(N_index,8); % KE + E1 + E2
    end            
    
    %% RMSE
    
    RMSE(j,1) = dt;
    RMSE(j,2) = sqrt( mse( sols(:,1) - refs(:,1) ) ); % q
    RMSE(j,3) = sqrt( mse( sols(:,2) - refs(:,2) ) ); % p
    RMSE(j,4) = sqrt( mse( sols(:,5) - refs(:,5) ) ); % S1 + S2
    RMSE(j,5) = sqrt( mse( sols(:,8) - E0 ) ); % KE + E1 + E2
     
    dt = dt*d_dt;
     
end
