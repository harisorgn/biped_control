% This script controls the search for limit cycles. 

clear

global d_sum % Accumulated step length, updated in biped3_step
global w_sum % Accumulated step width, updated in biped3_step
d_sum = 0 ;
w_sum = 0.0 ;

% Initialisation of the walking system as a structure
sys = struct('g',9.81,'k',0,'E',0,'v',0,'theta',0,'th0',0,'beta',0,'d',0,'w',0,'rot',0) ;

% Energy level and leg stiffness, their values would change during a
% bifurcation analysis.
sys.E = 1.0398 ;
sys.k = 17.8389 ;

% Initialisation of the vectors of limit cycle parameters. The column size
% is arbitrarily large, since the number of limit cycles is unknown at this
% point.
s_per = zeros(4,500) ; % phase space vector
eig_per = zeros(4,500) ; % eigenvalues of limit cycle
beta_per = zeros(1,500) ; % lateral (or azimuthal) leg angle relative to velocity vector
theta_per = zeros(1,500) ; % frontal (or polar) leg angle as measured from the horizontal plane
th0_per = zeros(1,500) ; % initial leg angle at TPO measured from the horizontal plane

i = 1 ;
% Range of initial leg angles at TPO
for th0 = 84 : 0.2 : 88
   
    sys.th0 = deg2rad(th0) ;
    
    % Solution of the 2x2 nonlinear system of the conditions for symmetric
    % walking
    q_init = [1, 0.32] ;
    options = optimoptions('fsolve','Display','off');
    f = @(q) biped3_init(q,sys) ;
    q_init = fsolve(f, q_init, options) ;
    
    % Initial phase space vector. Used as initial guess for Newton-Raphson.
    s0 = [real(q_init(1))*sin(sys.th0); real(q_init(1))*cos(sys.th0); 0; 0] ;
                
    % Range of theta angle
    for theta = 0.85 : 0.01 : 1.32
        
        sys.theta = theta ;
        
        % Range of beta angle
        for beta = 0.03 : 0.01 : 0.12
            
            sys.beta = beta ;
            
            % Call to Newton-Raphson function. The algorithm's numerical
            % parameters are set in the biped3_nr function.
            [s, eigenv, fall_chk] = biped3_nr(s0, sys) ;
            
            % Check for not falling and for the maximum absolute eigenvalue
            % being less than 1, as this is the condition for a stable
            % limit cycle. If the condition is true, then the limit cycle
            % is stable and its parameters are saved below.
            if (~fall_chk && max(abs(eigenv)) <= 1)
                
                s_per(:,i) = s ;
                eig_per(:,i) = eigenv ;
                th0_per(i) = th0 ;
                theta_per(i) = theta ;
                beta_per(i) = beta ;
                i = i + 1 ;
                
                fprintf('Limit cycle found for th0 = %f , beta = %f, theta = %f \n',[th0, beta, theta])
            end
        end
    end
end

% save('beta','s_per','eig_per','beta_per','theta_per','th0_per')
% exit


