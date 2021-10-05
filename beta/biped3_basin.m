% This function calculates two partial basins of attraction of a limit
% cycle s_lc, a 2D one in the y-z plane and a 1D one along perturbations of
% the polar velocity vector angle phi. The reason for splitting the calculation of the basin of
% attraction is described in the thesis report. However the calculation of 
% a 3D basin of attraction along the y, z and phi directions is encouraged.
% The 4th dimension of phase space, azimuthal velocity vector angle psi, was not included in the
% basin calculations, as both successful patterns were inherently neutrally
% stable to perturbations along that dimension.

function [s_basin, phi_range] = biped3_basin(s_lc, sys)

s_basin = zeros(1000, 4) ; % arbitrarily large row size
j = 1 ;

phi = s_lc(3) ;
psi = s_lc(4) ;

% Ranges of considered perturbations
y_min = 0.9 ;
y_max = 0.99 ;
z_min = 0.0 ;
z_max = 0.08 ;

n_steps = 50 ; % Number of necessary steps without fall

% Calculation of 2D basin
for y = y_min : 0.002 : y_max
    
    for z = z_min : 0.002 : z_max
        
        sys.v = sqrt(2*sys.E - sys.k*(1 - sqrt(y^2 + z^2))^2 - 2*y) ;
        q0 = [0; y; z; sys.v*cos(phi)*cos(psi); sys.v*sin(phi); sys.v*cos(phi)*sin(psi)] ;
        s = [q0(2); q0(3); asin(q0(5)/sqrt(q0(4)^2 + q0(5)^2 + q0(6)^2)); atan2(q0(6),q0(4))] ;
        
        fall_chk = 0 ;
        i = 1 ;
        
        while (~ fall_chk && i <= n_steps)
            
            [step, fall_chk] = biped3_step(q0, i, sys) ;
            q0 = step.q0' ;
            
            s = [q0(2); q0(3); asin(q0(5)/sqrt(q0(4)^2 + q0(5)^2 + q0(6)^2)); atan2(q0(6),q0(4))] ;
            
            i = i + 1 ;
            
        end
        
        % Euclidean norm is used to compare the values of y, z and phi of 
        % phase space between the 50th step and the limit cycle. 
        % Original similarity threshold set to 0.01
        if (~ fall_chk && norm(s(1:3) - s_lc(1:3)) < 0.01)
            
            s_basin(j,:) = [y, z, phi, psi] ;
            j = j + 1 ;
          
        end
        
    end
end

% Calculation of 1D basin
phi_range = zeros(200,1) ; % arbitrarily large row size
k = 1 ;

y = s_lc(1) ;
z = s_lc(2) ;

for phi = - 0.5 : 0.02 : 0.5
    
    sys.v = sqrt(2*sys.E - sys.k*(1 - sqrt(y^2 + z^2))^2 - 2*y) ;
    q0 = [sys.x_sect; y; z; sys.v*cos(phi)*cos(psi); sys.v*sin(phi); sys.v*cos(phi)*sin(psi)] ;
    s = [q0(2); q0(3); asin(q0(5)/sqrt(q0(4)^2 + q0(5)^2 + q0(6)^2)); atan2(q0(6),q0(4))] ;
    
    fall_chk = 0 ;
    i = 1 ;
    
    while (~ fall_chk && i <= n_steps)
        
        [step, fall_chk] = biped3_step(q0, i, sys) ;
        q0 = step.q0' ;
        
        s = [q0(2); q0(3); asin(q0(5)/sqrt(q0(4)^2 + q0(5)^2 + q0(6)^2)); atan2(q0(6),q0(4))] ;
        
        i = i + 1 ;
        
    end
    
    % The same similarity check is applied here as above
    if (~ fall_chk && norm(s(1:3) - s_lc(1:3)) < 0.01)
        
        phi_range(k) = phi ;
        k = k + 1 ;
      
    end
  
end

end