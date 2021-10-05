% Calculation of the double support dynamics during a step.

function dq = biped3_ds(t, q, sys)
   
    P = sys.k * (1 / sqrt(q(1)^2 + q(2)^2 + q(3)^2) - 1) ;
    Q = sys.k * (1 / sqrt((q(1) - sys.d)^2 + q(2)^2 + (q(3) - sys.w)^2) - 1) ; 
   
    dq = [ q(4) ;
           q(5) ;
           q(6) ;
           P * q(1) + Q * (q(1) - sys.d) ;
           P * q(2) + Q * q(2) - 1  ;
           P * q(3) + Q * (q(3) - sys.w)] ;
end