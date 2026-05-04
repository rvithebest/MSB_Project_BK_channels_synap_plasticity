clear;clc;close all;
% Parameter values and codes adapted from Aslam et al,2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.BASAL= 10;
params.BASAL1= 0.0;
params.BASAL2= 0.0;
params.k1= 0.000395*2.8;
params.k_1= 5*2.8;
params.k_2= 0.01*2.8;
params.k_4= 0.65* 0.0075*1.25*2.8;
params.k6= 0.982*2.8;
params.k_6= 0.1*2.8;
params.k7= 0.995*2.8;
params.k8 = 0.01*2.8;
params.k_8= 0.05*2.8;
params.k9= 0.09*2.8;
params.k10= 0.99*2.8;
params.k_10= 0.05*2.8;
params.k11= 0.85*2.8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params.kSYN1= 0.0058*2.8; % Default
params.kSYN1= 0.0058*2.11;  % Perturbation-1
% params.kSYN1= 0.0058*2.2;  % Perturbation-2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.k_12= 0.1*2.8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.lambda1= 0.00026*2.8; % Deafult
% params.lambda1= 0.00026*2.9;  % Perturbation-2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.P= 0.002214;
params.YT= 10.0;
params.T= 0.05;
params.k2= 0.008*2.8;
params.k4= 0.5*2.8;
params.k5= 0.022*2.8;
params.k3= 0.0075*2.8;
params.I2= 1.7495;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx=[1,2,3,4,5,6,7,8,9];
% idx=[1,2,3,4,5];
load('transformed_Ca_gatherer.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan = 0:1000:120000;
y00 = [10;0;0;0;0;0;0;0;0;0];
colors=jet(length(idx));
% colors- black, red, magenta, black, blue, violet
% colors=[0 0 0; 1 0 0; 1 0 1; 0 0 1; 0.5 0 0.5];
for i=1:length(idx)
    params.Ca=Ca_gatherer_transform{idx(i)};
    params.t_axis=t_gatherer{idx(i)};
    [t1, h1] = ode15s(@(t,y) eqn_solver(t,y,params), tspan, y00);
    xT1= h1(:,1)+ h1(:,2) + h1(:,3)+ 2*h1(:,5)+ 2*h1(:,6)+ h1(:,7) + h1(:,8);
    t11 = t1/60;
    F = h1(:,3)./xT1;
    subplot(1,2,1);
    hold on;
    plot(t11,xT1,'Color',colors(i,:),'LineWidth',2);
    xlabel('Time [min]');
    ylabel('Concentration of Total CaMKII [\muM]');
    xlim([0 2000]);
    ylim([0 30]);
    subplot(1,2,2);
    hold on;
    plot(t11,F,'Color',colors(i,:),'LineWidth',2);
    xlabel('Time [min]');
    ylabel('Fraction of CaMKII Active and Phosphorylated');
    xlim([0 2000]);
    ylim([0 0.6]);
end
% Legend: BK conductance values (BK_con_values=[0, 5e-4, 1e-3, 5e-3, 0.01, 0.015, 0.02, 0.025,0.03])
% legend('0','5e-4','1e-3','5e-3','0.01');
legend('0','5e-4','1e-3','5e-3','0.01','0.015','0.02','0.025','0.03');
