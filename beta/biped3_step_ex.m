clear

global d_sum
global w_sum
d_sum = 0 ;
w_sum = 0.0 ;

sys = struct('g',9.81,'k',0,'E',0,'v',0,'theta',0,'th0',0,'beta',0,'d',0,'w',0,'rot',0) ;

sys.E = 1.0398 ;
sys.k = 17.8389 ;

load('beta_lc.mat')

n_lc = 18 ;
n_step = 10 ;

q0 = q_init(n_lc,:) ;
sys.beta = beta(n_lc) ;
sys.theta = theta_abs(n_lc) ;

s_lc = [q0(2); q0(3); asin(q0(5)/sqrt(q0(4)^2 + q0(5)^2 + q0(6)^2)); atan2(q0(6),q0(4))] ;

q0 = q0 + [0 , 0 , 0 , 0 , 0 , 0] ;
q = [] ;
for i = 1 : n_step
    
    [step, fall_chk, ~] = biped3_step(q0, i, sys) ;
    
    if fall_chk
        break
    end
    
    q0 = step.q0 ;
    q = [q ; step.q] ;
        
end

figure(1)
plot(q(:,1),q(:,3),'b')
hold on
plot(q(:,1), q(:,3) + q(:,6)/sqrt(sys.k),'--g' ) % extrapolated centre of mass,
                                                 % assummed to be the spring leg equivalent of Hof's equation for stiff legs
xlabel('Forward position (m)')
ylabel('Lateral position (m)')
hold off

figure(2)
plot(q(:,1),q(:,2),'b')
xlabel('Forward position (m)')
ylabel('Vertical position (m)')

