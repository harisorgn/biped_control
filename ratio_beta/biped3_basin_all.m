clear

global d_sum
global w_sum
d_sum = 0 ;
w_sum = 0.0 ;

load('ratio_beta_lc.mat')

n = size(ratio,1) ;

ratio_beta_basin = zeros(500,2,n) ;
ratio_beta_phi = zeros(500,n) ;

for j = 1 : n
   
    sys = struct('g',9.81,'k',17.8389,'E',1.0398,'v',0,'theta',0,...
        'beta',beta(j),'ratio',ratio(j),'d',0,'w',0) ;
    
    q0 = q_init(j,:) ;
    
    s_lc = [q0(2); q0(3); asin(q0(5)/sqrt(q0(4)^2 + q0(5)^2 + q0(6)^2)); atan2(q0(6),q0(4))] ;
    
    [s_basin, phi_range] = biped3_basin(s_lc, sys) ;
    
    ratio_beta_basin(1:size(s_basin,1),:,j) = s_basin ;
    ratio_beta_phi(:,j) = phi_range ;
    
end

%save('ratio_beta_basin_all','ratio_beta_basin','ratio_beta_phi')
%exit