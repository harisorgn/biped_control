% Simulation of a single step, beginning from and ending at two consecutive
% instances of Transversal Plane Orientation (TPO, check thesis report).
% The variable i_step denotes the step number. Odd numbers denote left
% and even numbers right steps.

function [step, fall_chk, sys] = biped3_step(q0, i_step, sys)

    global d_sum
    global w_sum
    
    % Numerical parameters of integration algorithm
    rel_tol = 1e-12 ;
    abs_tol = 1e-12 ;
    max_step = 1e-3 ;

    % Step kinematics in a structure. 
    % ss1 is the first single support, ds is double support and ss2 the
    % second single support phase
    % q0 is the initial state vector and after the step it is the final
    % state vector or the initial vector of the next step. 
    % q is the accumulated state vector in global coordinates.
    step = struct('q',[],'q0',q0,'q_ss1',[],'q_ds',[],'q_ss2',[],'t_ss1',[],'t_ds',[],'t_ss2',[]) ;
    fall_chk = 0 ;
    
    t0 = 0 ;
    tf = 5 ;
    t_step = 0.01 ;
   
    opt=odeset('RelTol',rel_tol,'AbsTol',abs_tol,'MaxStep',max_step,'Events',@(t,q) biped3_TD(t,q,sys)) ;
    [step.t_ss1,q_ss,~,~,ie] = ode113(@(t,q) biped3_ss(t,q,sys), t0:t_step:tf, step.q0, opt) ;
    step.q0 = q_ss(end,:) ;
    
    if ie == 1 % check if single support was successful without violation of any of the falling conditions in biped3_TD
        
        % Calculation of strategy-specific leg angle parameters
        ang_v_lat = atan2(step.q0(6),step.q0(4)) ; % lateral (azimuthal) velocity angle
        ang_v_front = atan2(step.q0(5),sqrt(step.q0(4)^2 + step.q0(6)^2)) ; % frontal (polar) velocity angle
        
        sys.theta = sys.ratio*(pi/2 - abs(ang_v_front)) ; % theta angle relative to velocity vector
        sys.theta_abs = abs(-sys.theta + ang_v_front) ; % Absolute value of global frontal leg angle, as measured from the horizontal plane.
                                                        % Calculated in order to identify duplicate walking patterns in biped3_lc_param
                                                        
        sys.beta_abs = sys.beta + abs(ang_v_lat) ; % Absolute value of global lateral leg angle, as measured from the saggital plane
                                                   % Calculated in order to identify duplicate walking patterns in biped3_lc_param
        
        % Calculation of step length, width dependent on left/right step
        if mod(i_step,2) == 0
            sys.w = q_ss(end,3) + cos(-sys.theta + ang_v_front) * sin(-sys.beta + ang_v_lat) ;
            sys.d = q_ss(end,1) + cos(-sys.theta + ang_v_front) * cos(-sys.beta + ang_v_lat) ;
            
        else
            sys.w = q_ss(end,3) + cos(-sys.theta + ang_v_front) * sin(sys.beta + ang_v_lat) ;
            sys.d = q_ss(end,1) + cos(-sys.theta + ang_v_front) * cos(sys.beta + ang_v_lat) ;
            
        end
             
        opt=odeset('RelTol',rel_tol,'AbsTol',abs_tol,'MaxStep',max_step,'Events',@(t,q) biped3_LO(t,q,sys));
        [step.t_ds,q_ds,~,~,ie] = ode113(@(t,q) biped3_ds(t,q,sys), t0:t_step:tf, step.q0, opt) ;
        
        step.q0 = q_ds(end,:) - [sys.d 0 sys.w 0 0 0] ;
               
        step.q_ss1 = q_ss ;
        
        if ie ~= 2 % check if double support was successful without violation of any of the falling conditions in biped3_LO
            
            step.q = [step.q ; [d_sum 0 w_sum 0 0 0] + q_ss ; [d_sum 0 w_sum 0 0 0] + q_ds] ;
            step.q_ds = q_ds ;
            
            d_sum = d_sum + sys.d ; % update of global forward distance traveled
            w_sum = w_sum + sys.w ; % update of global lateral distance deviations. 
                                    % As step width alternates its sign during left/right steps, this value
                                    % is used to correctly plot travel paths, when step widths differ
                                    % between left and right legs
            
            opt=odeset('RelTol',rel_tol,'AbsTol',abs_tol,'MaxStep',max_step,'Events',@(t,q) biped3_TPO(t,q, sys));
            [step.t_ss2,q_ss,~,~,ie] = ode113(@(t,q) biped3_ss(t,q,sys), t0:t_step:tf, step.q0, opt) ;
            step.q0 = q_ss(end,:) ;
            
            if isempty(ie)
                ie = 0 ;
            end
            
            if (ie == 2 ) % check if second single support was successful without violation of any of the falling conditions in biped3_TPO
                fall_chk = 1 ;
            else
                step.q = [step.q ; [d_sum 0 w_sum 0 0 0] + q_ss] ;
                step.q_ss2 = q_ss ;
            end
        else
            fall_chk = 1 ;
        end
    else
        fall_chk = 1 ;
    end
end