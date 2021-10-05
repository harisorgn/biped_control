% Simulation of one stride. This function is used by the biped3_nr, where
% 4 strides have to be performed to calculate all partial derivatives of
% the Poincare map Jacobian. It was developed to make biped3_nr more
% concise.

function [stride, fall_chk] = biped3_stride(q0, sys)

    stride = struct('q',[],'q0',[]) ;
    
    [step1, fall_chk, ~] = biped3_step(q0, 1, sys) ;
    
    if ~ fall_chk
       
        [step2, fall_chk, ~] = biped3_step(step1.q0, 2, sys) ;
        
        if ~ fall_chk
           
            stride.q = [step1.q ; step2.q] ;
            stride.q0 = step2.q0 ;
            
        end
        
    end
end