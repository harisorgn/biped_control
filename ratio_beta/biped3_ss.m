% Calculation of the single support dynamics during a step.

function dq = biped3_ss(t, q, sys)
    
    P = sys.k * (1 / sqrt(q(1)^2 + q(2)^2 + q(3)^2) - 1) ;
    
    dq = [ q(4) ;
           q(5) ;
           q(6) ;
           P * q(1) ;
           P * q(2) - 1 ;
           P * q(3) ] ;
end