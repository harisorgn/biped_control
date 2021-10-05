clear

global d_sum
global w_sum
d_sum = 0 ;
w_sum = 0.0 ;

sys = struct('g',9.81,'k',0,'E',0,'v',0,'theta',0,'th0',0,'beta',0,'d',0,'w',0,'ratio',0,'rot',0) ;

sys.E = 1.0398 ;
sys.k = 17.8389 ;

s_per = zeros(4,500) ;
eig_per = zeros(4,500) ;
beta_per = zeros(1,500) ;
theta_per = zeros(1,500) ;
th0_per = zeros(1,500) ;
ratio_per = zeros(1,500) ;


load('ratio_beta_lc.mat')

n = 11 ;

q0 = q_init(n,:) ;
sys.beta = beta(n) ;
sys.ratio = ratio(n) ;

s_lc = [q0(2); q0(3); asin(q0(5)/sqrt(q0(4)^2 + q0(5)^2 + q0(6)^2)); atan2(q0(6),q0(4))] ;

[ratio_beta_basin, ratio_beta_phi] = biped3_basin(s_lc, sys) ;

% save('ratio_beta_basin_ex','ratio_beta_basin','ratio_beta_phi')
% exit