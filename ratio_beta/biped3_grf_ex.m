clear

global d_sum
global w_sum
d_sum = 0 ;
w_sum = 0.0 ;

sys = struct('g',9.81,'k',0,'E',0,'v',0,'theta',0,'th0',0,'beta',0,'d',0,'w',0,'ratio',0,'rot',0) ;

sys.E = 1.0398 ;
sys.k = 17.8389 ;

load('ratio_beta_lc.mat')

n = 10 ;

q0 = q_init(n,:) ;
sys.beta = beta(n) ;
sys.ratio = ratio(n) ;

[fr, fl] = biped3_grf(q0,sys) ;