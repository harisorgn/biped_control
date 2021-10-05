% Calculation of the ground reaction force profiles along x,y and z for the left and right
% legs, fl and fr respectively.
% The function plots the two force profiles along the vertical direction y
% and returns two Nx3 matrices, where each column represents the force
% profile of one of the 3 spatial directions (x,y,z) during each of the N time
% steps.

% Forces are present as long as the leg is in contact with the ground, so
% each force calculation begins during the initiation of double support and
% ends when the the next double support phase ends.

function [fr, fl] = biped3_grf(q0,sys)

    [step1, ~, sys] = biped3_step(q0, 1, sys) ;
    d1 = sys.d ;
    w1 = sys.w ;
    [step2, ~, sys] = biped3_step(step1.q0, 2, sys) ;
    d2 = sys.d ;
    w2 = sys.w ;
    [step3, ~, sys] = biped3_step(step2.q0, 3, sys) ;

    fr = [] ;
    fl = [] ;
    
    for i = 1 : size(step1.q_ds,1)
       
        Q = sys.k * (1 / sqrt((step1.q_ds(i,1) - d1)^2 + step1.q_ds(i,2)^2 + (step1.q_ds(i,3) - w1)^2) - 1) ; 
        
        fr = [fr ; Q * (step1.q_ds(i,1) - d1), Q * step1.q_ds(i,2), Q * (step1.q_ds(i,3) - w1)] ;
        
    end
    
    for i = 2 : size(step1.q_ss2,1)
       
        P = sys.k * (1 / sqrt(step1.q_ss2(i,1)^2 + step1.q_ss2(i,2)^2 + step1.q_ss2(i,3)^2) - 1) ;
        
        fr = [fr ; P * step1.q_ss2(i,1), P * step1.q_ss2(i,2), P * step1.q_ss2(i,3)] ;
        
    end
    
    for i = 2 : size(step2.q_ss1,1) - 1
       
        P = sys.k * (1 / sqrt(step2.q_ss1(i,1)^2 + step2.q_ss1(i,2)^2 + step2.q_ss1(i,3)^2) - 1) ;
        
        fr = [fr ; P * step2.q_ss1(i,1), P * step2.q_ss1(i,2), P * step2.q_ss1(i,3)] ;
        
    end
        
    for i = 1 : size(step2.q_ds,1)
       
        P = sys.k * (1 / sqrt(step2.q_ds(i,1)^2 + step2.q_ds(i,2)^2 + step2.q_ds(i,3)^2) - 1) ;
        
        fr = [fr ; P * step2.q_ds(i,1), P * step2.q_ds(i,2), P * step2.q_ds(i,3)] ;
        
        Q = sys.k * (1 / sqrt((step2.q_ds(i,1) - d2)^2 + step2.q_ds(i,2)^2 + (step2.q_ds(i,3) - w2)^2) - 1) ; 
        
        fl = [fl ; Q * (step2.q_ds(i,1) - d2), Q * step2.q_ds(i,2), Q * (step2.q_ds(i,3) - w2)] ;
    end
    
    for i = 2 : size(step2.q_ss2,1)
       
        P = sys.k * (1 / sqrt(step2.q_ss2(i,1)^2 + step2.q_ss2(i,2)^2 + step2.q_ss2(i,3)^2) - 1) ;
        
        fl = [fl ; P * step2.q_ss2(i,1), P * step2.q_ss2(i,2), P * step2.q_ss2(i,3)] ;
        
    end
    
    for i = 2 : size(step3.q_ss1,1)
       
        P = sys.k * (1 / sqrt(step3.q_ss1(i,1)^2 + step3.q_ss1(i,2)^2 + step3.q_ss1(i,3)^2) - 1) ;
        
        fl = [fl ; P * step3.q_ss1(i,1), P * step3.q_ss1(i,2), P * step3.q_ss1(i,3)] ;
        
    end
    
    for i = 2 : size(step3.q_ds,1)
       
        P = sys.k * (1 / sqrt(step3.q_ds(i,1)^2 + step3.q_ds(i,2)^2 + step3.q_ds(i,3)^2) - 1) ;
        
        fl = [fl ; P * step3.q_ds(i,1), P * step3.q_ds(i,2), P * step3.q_ds(i,3)] ;
        
    end
    
    x_axr = linspace(0,1,size(fr,1)) ;
    x_axl = linspace(0,1,size(fl,1)) ;
    
    subplot(1,2,1)
    plot(x_axr,fr(:,2))
    title('Right leg')
    ylabel('Ground reaction force (body weight)')
    
    set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',16,'FontName','Arial')
    set(findall(gcf,'type','text'),'FontUnits','points','FontWeight','normal','FontSize',16,'FontName','Arial')
    
    subplot(1,2,2)
    plot(x_axl,fl(:,2))
    title('Left leg')
    ylabel('Ground reaction force (body weight)')
    
    set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',16,'FontName','Arial')
    set(findall(gcf,'type','text'),'FontUnits','points','FontWeight','normal','FontSize',16,'FontName','Arial')
end