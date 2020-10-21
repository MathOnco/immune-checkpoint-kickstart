clc;clear;close all;
close all;
lw1 = 3;
colors = get(groot,'DefaultAxesColorOrder');
options = odeset('RelTol',1e-9,'AbsTol',1e-9);


% parameters
eps = 0.01;
lambda = 0.5;
kT = 0.4;
b2 = 0.1;
c1 = 0.1; % CCR7 cost
c2 = 0.1; % PD-L1 cost
k1 = 1.1;
k2 = 0.95;
ff1 = 0.2; % initial condition of PD-L1
ff2 = 0.1; % initial condition of CCR7
v0 = 1;
Q = eye(4); % no mutations

totalT = 6; % 6 months
x0 = [(1 - ff1)*(1 - ff2), ff1*(1-ff2), ff2*(1-ff1), ff1*ff2, ...
           v0*(1 - ff1)*(1 - ff2), v0*ff1*(1-ff2), v0*ff2*(1-ff1), v0*ff1*ff2];

% each row is a new scenario 
%  --> (col 1 = AI; col 2 = anti-PDL1)
treatments = [0,0;
              1,0;
              0,1;
              1,1];
       
for scenario = 1:4
    
    drug1 = treatments(scenario,1);
    drug2 = treatments(scenario,2);
    figure(scenario);
        
    A = payoff(kT, k1*drug1, k2*drug2, b2, c1, c2);
    K1 = k1*drug1;
    K2 = k2*drug2;


    [t, xx]=ode45(@(t,n)rep_ode(t, n, A, Q, lambda), [0 totalT], x0,options);

    % determine the winner:
    xF = xx(end,1:4);
    [~,colorUnderi] = max(xF);

    % fill color underneath the winner:
    a = patch([t' flip(t')], [xx(:,colorUnderi)' (zeros(1,length(t)))],'k'); hold on;
    a.FaceColor = colors(colorUnderi,:);
    a.FaceAlpha = 0.2;
    a.EdgeColor = 'none';

    % plot all colors:
    plot(t, xx(:,1),'-','LineWidth',lw1,'Color',colors(1,:)); hold on;
    plot(t, xx(:,2),'-','LineWidth',lw1,'Color',colors(2,:)); hold on;
    plot(t, xx(:,3),'-','LineWidth',lw1,'Color',colors(3,:)); hold on;
    plot(t, xx(:,4),'-','LineWidth',lw1,'Color',colors(4,:)); hold on;

    xlim([0 totalT]);
    ylim([0 1]);

    %% clean replicator dynamics figures
    nice_plot(scenario,'time (months)','$x_i(t)$',true);

end

