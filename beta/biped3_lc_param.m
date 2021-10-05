% This script identifies unique limit cycles among all the limit cycles
% found by biped3_nr through biped3_lc_search and calculates useful parameters of their gaits.

% Since both successful strategies, ratio-beta and beta, are neutrally
% stable to lateral perturbations, the same limit cycle pattern may be found
% more than once by the Newton-Raphson method having walking directions at different angles 
% with respect to the global forward x axis. 

% The walking direction angles are found and limit cycle patterns are
% rotated to align the walking direction with the global x axis in order to
% identify and duplicates.


clear
load('beta.mat')

global d_sum
global w_sum
d_sum = 0 ;
w_sum = 0.3 ;

n = 1 ;
while (th0_per(n) ~= 0)
    n = n + 1 ;
end

dl = zeros(n-1,1) ; % left step length
wl = zeros(n-1,1) ; % left step width

dr = zeros(n-1,1) ; % right step length
wr = zeros(n-1,1) ; % right step width

t_sl = zeros(n-1,1) ; % left step time
t_sr = zeros(n-1,1) ; % right step time

v_fwd = zeros(n-1,1) ; % mean forward velocity
v_lat = zeros(n-1,1) ; % mean lateral velocity
v_lat_TD = zeros(n-1,1) ; % lateral velocity at touch-down event
v_fwd_TD = zeros(n-1,1) ; % forward velocity at touch-down event

q_init = zeros(n-1, 6) ;
beta = zeros(n-1,1) ;
theta_abs = zeros(n-1,1) ;
beta_abs = zeros(n-1,1) ;

r = 1 ;

for j = 1 : n-1 
    
    sys = struct('g',9.81,'k',17.8389,'E',1.0398,'v',0,'theta',theta_per(j),'th0',th0_per(j),'beta',beta_per(j),...
                'd',0,'w',0,'theta_abs',0,'beta_abs',0) ;

    s0 = s_per(:,j) ;
    sys.v = sqrt(2*sys.E - sys.k*(1 - sqrt(s0(1)^2 + s0(2)^2))^2 - 2*s0(1)) ;
    
    q0 = [0; s0(1); s0(2); sys.v*cos(s0(3))*cos(s0(4)); sys.v*sin(s0(3)); sys.v*cos(s0(3))*sin(s0(4))] ;
    
    % One stride to find the walking direction and rotate the reference frame accordingly
    q = [] ;
    for i = 1 : 2
        
        [step, fall_chk, ~] = biped3_step(q0, i, sys) ;
        
        q0 = step.q0 ;
        q = [q ; step.q] ;
    end
    
    sys.rot = - atan2(q(end,3) - q(1,3), q(end,1) - q(1,1)) ;
    
    q0 = [-s0(2)*sin(sys.rot) ; s0(1) ; s0(2)*cos(sys.rot) ; sys.v*cos(s0(3))*cos(s0(4) + sys.rot) ; sys.v*sin(s0(3)) ; sys.v*cos(s0(3))*sin(s0(4) + sys.rot)] ;
    
    % One stride to find the state q_init at TPO (x = 0) to be used as the 
    % initial state q0 for integration, now that the walking direction coincides 
    % with the global x axis
    for i = 1 : 2
        
        [step, fall_chk, ~] = biped3_step(q0, i, sys) ;
        
        q0 = step.q0 ;
    end
    
    % Three steps to calculate left/right step length, width, time, velocities only for unique patterns
    [step1, ~, sys] = biped3_step(q0, 1, sys) ;
    
    % Check for duplicates
    pattern_in = 0 ;
    for k = 1 : j
       if (abs(theta_abs(k) - sys.theta) < 1e-5 && abs(beta_abs(k) - sys.beta_abs) < 1e-5 ...
           && abs(beta(k) - sys.beta) < 1e-5 )
       
           pattern_in = 1 ;
           break
       end
    end
    
    % Calculation of gait parameters for unique patterns
    if ~ pattern_in
        q_init(r, :) = q0 ;
        
        dr(r) = sys.d ;
        wr(r) = sys.w ;
    
        beta(r) = sys.beta ;
        theta_abs(r) = sys.theta ;
        beta_abs(r) = sys.beta_abs ;
        
        q0 = step1.q0 ;
        [step2, ~, sys] = biped3_step(q0, 2, sys) ;
        
        dl(r) = sys.d ;
        wl(r) = sys.w ;
        
        t_sr(r) = step1.t_ds(end) - step1.t_ds(1) + step1.t_ss2(end) - step1.t_ss2(1) + step2.t_ss1(end) - step2.t_ss1(1) + step2.t_ds(end) - step2.t_ds(1) ;
        
        q0 = step2.q0 ;
        [step3, ~, sys] = biped3_step(q0, 3, sys) ;
        
        t_sl(r) = step2.t_ds(end) - step2.t_ds(1) + step2.t_ss2(end) - step2.t_ss2(1) + step3.t_ss1(end) - step3.t_ss1(1) + step3.t_ds(end) - step3.t_ds(1) ;
        
        v_fwd(r) = mean([step1.q(:,4) ; step2.q(:,4)]) ;
        v_lat(r) = mean([abs(step1.q(:,6)) ; abs(step2.q(:,6))]) ;
        v_lat_TD(r) = mean([abs(step1.q_ds(1,6)) ; abs(step2.q_ds(1,6))]) ;
        v_fwd_TD(r) = mean([step1.q_ds(1,4) ; step2.q_ds(1,4)]) ;
        
        r = r+1 ;
    end  
end

n = 1 ;
while (t_sl(n) ~= 0)
    n = n + 1 ;
end

dl = dl(1:n-1) ;
wl = wl(1:n-1) ;

dr = dr(1:n-1) ;
wr = wr(1:n-1) ;

t_sl = t_sl(1:n-1) ;
t_sr = t_sr(1:n-1) ;

v_fwd = v_fwd(1:n-1) ;
v_lat = v_lat(1:n-1) ;
v_lat_TD = v_lat_TD(1:n-1) ;
v_fwd_TD = v_fwd_TD(1:n-1) ;

q_init = q_init(1:n-1,:) ;
beta = beta(1:n-1) ;
theta_abs = theta_abs(1:n-1) ;
beta_abs = beta_abs(1:n-1) ;

save('beta_lc','q_init','dl','wl','dr','wr','t_sl','t_sr','beta','theta_abs','beta_abs','v_fwd','v_lat','v_lat_TD','v_fwd_TD')
