% The lift-off event

function [value,isterminal,direction] = biped3_LO(~, q, ~)
    
    value = [sqrt(q(1)^2 + q(2)^2 + q(3)^2) - 1 ; q(2) ; q(4)] ;
    isterminal = [1 ; 1 ; 1] ;
    direction = [1 ; -1 ; -1] ;
    
end