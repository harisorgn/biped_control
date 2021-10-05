% Calculation of the initial leg length q(1) and velocity magnitude q(2),
% by solving a 2x2 nonlinear system.

% This function is called by biped3_lc_search and is solved by the built-in
% Matlab fsolve function, in order to set initial conditions for symmetric
% walking.

function F = biped3_init(q, sys)

    F(1) = 1 - q(1) - (sys.k + sin(sys.th0) - sqrt(sys.k^2 + 2*sys.k*sin(sys.th0) + sin(sys.th0)^2 - 4*sys.k*sin(sys.th0) + 4*sys.k*q(2)^2))/(2*sys.k) ;
    F(2) = q(2) - sqrt(2*sys.E - sys.k*(1 - q(1))^2 - 2*q(1)*sin(sys.th0)) ;
    
end