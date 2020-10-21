clc;clear;close all;
close all;
ms = 60;
lw = 1;

options = odeset('RelTol',1e-9,'AbsTol',1e-9);

% parameters
lambda = 0.5;
kT = 0.4;
b2 = 0.1;
c1 = 0.1;
c2 = 0.1;
k1 = 1.1;
k2 = 0.95;
ff1 = 0.1; % initial condition of PD-L1
ff2 = 0.1; % initial condition of CCR7
v0 = 1;
totalT = 6; % months
Q = eye(4);

%% generate color gradient
figure(1);hold on;
scale = 100;

% metastasis background gradient:
p1 = 0.5; p2 = 0.5;
c = gradientOfMetPotential(p1,p2);
h = pcolor(c); hold on;
shading interp;
set(h, 'EdgeColor', 'none');
clean_plot(1,'','',true);
colorMap = [linspace(1,0,256)', linspace(1,0,256)',linspace(1,0,256)'];%zeros(256,2)];
colormap(colorMap);
colorbar; axis ij; axis square;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
xlim([1, scale]);
ylim([1, scale]);
box on;
set(gca,'linewidth',6);

% line colors:
colors = [0,0,0; 0.0784,0.0039,0.8000; 1,0,0;0.5547,0.3984,0.7266];


% each row is a new scenario 
%  --> (col 1 = AI; col 2 = anti-PDL1)
treatments = [0,0;
              1,0;
              0,1;
              1,1];
       
for scenario = 1:4
    
    drug1 = treatments(scenario,1);
    drug2 = treatments(scenario,2);
    x0 = [(1 - ff1)*(1 - ff2), ff1*(1-ff2), ff2*(1-ff1), ff1*ff2, v0*(1 - ff1)*(1 - ff2), v0*ff1*(1-ff2), v0*ff2*(1-ff1), v0*ff1*ff2];

    A = payoff(kT, k1*drug1, k2*drug2, b2, c1, c2);
    K1 = k1*drug1;
    K2 = k2*drug2;
    [tt, xx]=ode45(@(t,n)rep_ode(t, n, A, Q, lambda), [0 totalT], x0,options);
    CCR7 = (xx(:,3) + xx(:,4));
    PDL1 = (xx(:,2) + xx(:,4));


    %% figure 3B
    figure(1); hold on;
    plot(CCR7*scale + 1,PDL1*scale + 1,'-', 'LineWidth', lw,'Color', colors(scenario,:)); hold on;
    plot(CCR7(end)*scale + 1,PDL1(end)*scale + 1,'.','MarkerSize',ms,'Color',colors(scenario,:)); hold on;
    
    % add the black starting dot
    x0 = [(1 - ff1)*(1 - ff2), ff1*(1-ff2), ff2*(1-ff1), ff1*ff2];
    CCR70 = (x0(3) + x0(4));
    PDL10 = (x0(2) + x0(4));
    plot(CCR70*scale + 1,PDL10*scale + 1,'.', 'MarkerSize', ms,'Color', 'Black'); hold on;


end











