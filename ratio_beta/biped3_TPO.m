% The Transversal Plane Orientation event. This is the initial and the
% final state during each step.

function [value,isterminal,direction] = biped3_TPO(~, q, ~)
    
    value = [q(1) ; q(2) ; q(4)] ;
    
    isterminal = [1 ; 1 ; 1] ;
    
    direction = [1 ; -1 ; -1] ;
    
end