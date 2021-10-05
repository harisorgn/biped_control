% Plotting function for the calculated 2D basins of attraction along y and z 
% with respect to the corresponding control variables (ie leg angles theta and beta). 
% It works both for plotting the basin of a single pattern, as is the output of
% biped3_basin_ex and for plotting the basins of all patterns of a
% strategy, as provided by the biped3_basin_all function.

clear
load('ratio_beta_basin_all.mat') % Change this to ratio_beta_basin_ex.mat for plotting the basin of a single pattern
load('ratio_beta_lc.mat')

count = zeros(size(ratio_beta_basin,3),1) ;
bound = zeros(30, size(ratio_beta_basin,3)) ;
bound_sz = zeros(size(ratio_beta_basin,3),1) ;

for k = 1:size(ratio_beta_basin,3)
    for i = 1:size(ratio_beta_basin,1)
        
        if ratio_beta_basin(i,1,k) ~= 0
            count(k) = count(k) + 1 ;
        end
        
    end
    k_bound = boundary(ratio_beta_basin(1:count(k),:,k),0) ;
    bound_sz(k) = size(k_bound,1) ;
    bound(1:bound_sz(k),k) = k_bound ;
end

% Colour scale bounds
ran=range(theta_abs) ; 
min_val=min(theta_abs) ;
max_val=max(theta_abs) ;
col_val = floor(((theta_abs - min_val)/ran)*63)+1; 
col = zeros(size(ratio_beta_basin,3),3) ;
area = zeros(size(ratio_beta_basin,3),1) ;

% Colour map choice
p = colormap('winter') ;

line_sn = [beta_abs, q_init(:,3), q_init(:,2)] ;

figure(1)
for k = 1:size(ratio_beta_basin,3)
   
    col(k,:) = p(col_val(k),:) ;
    
    y = ratio_beta_basin(bound(1:bound_sz(k),k),1,k) ; 
    z = ratio_beta_basin(bound(1:bound_sz(k),k),2,k) ; 
    theta_plot = ones(bound_sz(k),1) * theta_abs(k) ;
    beta_plot = ones(bound_sz(k),1) * beta_abs(k) ;
    
    area(k) = polyarea(z,y) ;
    
    plot3(beta_plot,z,y,'color',col(k,:))
    hold on
    plot3(line_sn(k,1), line_sn(k,2), line_sn(k,3),'*','color',col(k,:))
end

caxis([min_val max_val])
xlabel('Lateral leg angle \beta (rad)')
ylabel('Lateral position z')
zlabel('Vertical Position y')
%title('Basin of attraction of patterns of ratio-beta strategy')
grid on

cb = colorbar ;
ylabel(cb,'Frontal leg angle \theta (rad)','FontWeight','normal','FontSize',16,'FontName','Arial')
hold off

set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',16,'FontName','Arial')
set(findall(gcf,'type','text'),'FontUnits','points','FontWeight','normal','FontSize',16,'FontName','Arial')

if (size(ratio_beta_basin,3) > 1)
    
    figure(2)
    ran=range(area) ;
    min_val=min(area) ;
    max_val=max(area) ;
    col_val = floor(((area - min_val)/ran)*63)+1;
    
    p = colormap('winter') ;
    col = p(col_val,:) ;
    
    for k = 1:size(ratio_beta_basin,3)
        
        plot(beta_abs(k),theta_abs(k),'o','color',col(k,:))
        hold on
    end
    
    caxis([min_val max_val])
    xlabel('Lateral leg angle \beta (rad)')
    ylabel('Frontal leg angle \theta (rad)')
    %title('Area of basin of attraction for each pattern')
    grid on
    
    cb = colorbar ;
    ylabel(cb,'Area of basin of attraction','FontWeight','normal','FontSize',16,'FontName','Arial')
    hold off
    
    set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',16,'FontName','Arial')
    set(findall(gcf,'type','text'),'FontUnits','points','FontWeight','normal','FontSize',16,'FontName','Arial')

end