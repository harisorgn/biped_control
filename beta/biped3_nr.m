% Implementation of the Newton-Raphson root-finding method.

function [s_root, eigenv, fall_chk] = biped3_nr(s0, sys)

    fall_chk = 0 ; % fall check
    eigen_calculated = 0 ; % after the root is found, one more loop is performed 
                           % to calculate the limit cycle's Jacobian and
                           % its eigenvalue. This variables checks if this
                           % last loop has been completed before exiting
                           % the function.
    eigen_chk = 0 ; % This variable is used to switch between the normal Newton-Raphson 
                    % iterations and the calculation of the limit cycle
                    % Jacobian during the last loop.
    
    eps = 1e-6 ; % convergence accuracy
    fd_eps = 1e-8 ; % step size for the approximation of derivatives by finite forward differences in the Jacobian

    s_root = zeros(4,1) ; % the root of the x(i+1) = P(x(i)) equation corresponding to a limit cycle, where P is the Poincare map
    eigenv = zeros(4,1) ; % eigenvalues of the Poincare map Jacobian
    
    s_next = s0 ;
    
    i = 1 ; % iteration counter
    max_iter = 200 ; % maximum number of iterations, set as a fail-safe to the Newton-Raphson algorithm
    
    while (~fall_chk && ~eigen_calculated && i <= max_iter) % Loops run until the eigenvalues 
                                                            % of the Poincare map Jacobian are calculated, unless the system
                                                            % falls before that.
        
        s_current = s_next ;
        
        v_current = sqrt(2*sys.E - sys.k*(1 - sqrt(s_current(1)^2 + s_current(2)^2))^2 - 2*s_current(1)) ;
        
        q_current = [0, s_current(1), s_current(2), ... 
            v_current*cos(s_current(3))*cos(s_current(4)), v_current*sin(s_current(3)), v_current*cos(s_current(3))*sin(s_current(4))] ;
                
        [stride, fall_chk] = biped3_stride(q_current, sys) ;
        q_map = stride.q0 ;
        
        if ~ fall_chk
            
            s_map = [q_map(2); q_map(3); asin(q_map(5)/sqrt(q_map(4)^2 + q_map(5)^2 + q_map(6)^2)); atan2(q_map(6),q_map(4))] ;
        
            
            v_dy = sqrt(2*sys.E - sys.k*(1 - sqrt((s_current(1) + fd_eps)^2 + s_current(2)^2))^2 - 2*(s_current(1) + fd_eps)) ;
            
            q_dy = [0, s_current(1) + fd_eps, s_current(2), ...
                v_dy*cos(s_current(3))*cos(s_current(4)), v_dy*sin(s_current(3)), v_dy*cos(s_current(3))*sin(s_current(4))] ;
            
            [stride, fall_chk] = biped3_stride(q_dy, sys) ;
            
        if ~ fall_chk
            q_dy = stride.q0 ;
            s_dy = [q_dy(2); q_dy(3); asin(q_dy(5)/sqrt(q_dy(4)^2 + q_dy(5)^2 + q_dy(6)^2)); atan2(q_dy(6),q_dy(4))] ;
            
            v_dz = sqrt(2*sys.E - sys.k*(1 - sqrt(s_current(1)^2 + (s_current(2) + fd_eps)^2))^2 - 2*s_current(1)) ;
            
            q_dz = [0, s_current(1), s_current(2) + fd_eps, ...
                v_dz*cos(s_current(3))*cos(s_current(4)), v_dz*sin(s_current(3)), v_dz*cos(s_current(3))*sin(s_current(4))] ;
            
            [stride, fall_chk] = biped3_stride(q_dz, sys) ;

        if ~ fall_chk
            q_dz = stride.q0 ;
            s_dz = [q_dz(2); q_dz(3); asin(q_dz(5)/sqrt(q_dz(4)^2 + q_dz(5)^2 + q_dz(6)^2)); atan2(q_dz(6),q_dz(4))] ;
             
            q_dphi = [0, s_current(1), s_current(2), ...
                v_current*cos(s_current(3) + fd_eps)*cos(s_current(4)), v_current*sin(s_current(3) + fd_eps), v_current*cos(s_current(3) + fd_eps)*sin(s_current(4))] ;
            
            [stride, fall_chk] = biped3_stride(q_dphi, sys) ;
            
        if ~ fall_chk
            q_dphi = stride.q0 ;
            s_dphi = [q_dphi(2); q_dphi(3); asin(q_dphi(5)/sqrt(q_dphi(4)^2 + q_dphi(5)^2 + q_dphi(6)^2)); atan2(q_dphi(6),q_dphi(4))] ;
            
            q_dpsi = [0, s_current(1), s_current(2), ...
                v_current*cos(s_current(3))*cos(s_current(4) + fd_eps), v_current*sin(s_current(3)), v_current*cos(s_current(3))*sin(s_current(4) + fd_eps)] ;
            
            [stride, fall_chk] = biped3_stride(q_dpsi, sys) ;
            
        if ~fall_chk
            q_dpsi = stride.q0 ;
            s_dpsi = [q_dpsi(2); q_dpsi(3); asin(q_dpsi(5)/sqrt(q_dpsi(4)^2 + q_dpsi(5)^2 + q_dpsi(6)^2)); atan2(q_dpsi(6),q_dpsi(4))] ;
            
            F = s_current - s_map ; % F = 0 represents the condition for limit cycles
            
            % Jacobian of the Poincare map
            J = [(s_dy(1) - s_map(1))/fd_eps , (s_dz(1) - s_map(1))/fd_eps , (s_dphi(1) - s_map(1))/fd_eps , (s_dpsi(1) - s_map(1))/fd_eps ; ...
                (s_dy(2) - s_map(2))/fd_eps , (s_dz(2) - s_map(2))/fd_eps , (s_dphi(2) - s_map(2))/fd_eps , (s_dpsi(2) - s_map(2))/fd_eps ; ...
                (s_dy(3) - s_map(3))/fd_eps , (s_dz(3) - s_map(3))/fd_eps , (s_dphi(3) - s_map(3))/fd_eps , (s_dpsi(3) - s_map(3))/fd_eps ; ...
                (s_dy(4) - s_map(4))/fd_eps , (s_dz(4) - s_map(4))/fd_eps , (s_dphi(4) - s_map(4))/fd_eps , (s_dpsi(4) - s_map(4))/fd_eps ] ;
       
            % Calculation of the next phase space vector by Newton-Raphson
            s_next = s_current - (eye(4) - J)\F ; 
            
            % Euclidean norm is used to compare two consecutive root approximations
            converge_chk = norm(s_next - s_current) ; 
          
            % The velocity magnitude of the next iteration's phase space
            % vector is calculated from the energy conservation equation.
            v_next = sqrt(2*sys.E - sys.k*(1 - sqrt(s_next(1)^2 + s_next(2)^2))^2 - 2*s_next(1)) ;
            
            % The check for falling includes invalid values for the next
            % iteration's phase space vector.
            if (sqrt(s_next(1)^2 + s_next(2)^2) > 1 || s_next(1) < 0 || s_next(2) < 0 || cos(s_next(3))*cos(s_next(4)) < 0 || ~isreal(v_next))
               
                fall_chk = 1 ;
                
            end
            
            if (converge_chk <= eps || eigen_chk > 0)
                
                eigen_chk = eigen_chk + 1 ;
                if (eigen_chk == 1)
                    % The first time this block is accessed is when the method
                    % has converged to a root. The corresponding phase space vector is
                    % saved and the next iteration vector remains equal to
                    % the root.
                    s_root = s_next ;
                    s_next = s_root + 0.00 * s_root ;
                elseif (eigen_chk == 2)
                    % The second time this block is accessed is after the
                    % calculation of the Jacobian by finite differences
                    % while having the found root as the initial vector.
                    eigen_calculated = 1 ;
                    eigenv = eig(J) ;
                end
                
            end          
        end
        end
        end
        end
        end
    end
    

end
