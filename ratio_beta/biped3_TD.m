% The touch-down event. This is a strategy-specific function.

function [value,isterminal,direction] = biped3_TD(~, q, sys)
    
    ang_v_front = atan2(q(5),sqrt(q(4)^2 + q(6)^2)) ; % frontal (polar) velocity angle
    
    sys.theta = sys.ratio*(pi/2 - abs(ang_v_front)) ; % theta angle relative to velocity vector
    
    value = [q(2) - sin(sys.theta + abs(ang_v_front)) ; q(2) ; q(4) ; sqrt(q(1)^2 + q(2)^2 + q(3)^2) - 1] ;
    
    isterminal = [1 ; 1 ; 1 ; 1] ;
    direction = [-1 ; -1 ; -1 ; 1] ;

end