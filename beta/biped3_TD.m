% The touch-down event. This is a strategy-specific function.

function [value,isterminal,direction] = biped3_TD(~, q, sys)
            
    value = [q(2) - sin(sys.theta) ; q(2) ; q(4) ; sqrt(q(1)^2 + q(2)^2 + q(3)^2) - 1] ;
    
    isterminal = [1 ; 1 ; 1 ; 1] ;
    direction = [-1 ; -1 ; -1 ; 1] ;

end