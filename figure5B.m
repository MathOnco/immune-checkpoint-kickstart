clc;clear;close all;
close all;

mss = 45;
lw = 1;
colors = get(groot,'DefaultAxesColorOrder');    
options = odeset('RelTol',1e-9,'AbsTol',1e-9);


%% delayed AI 6 months vs delayed anti-PDL1 6 months
purple = [0.5547    0.3984    0.7266];
red = [1,0,0];
blue = [0.0784    0.0039    0.8000];
colors1 = ColorRange(blue,purple,7);
colors2 = ColorRange(purple,red,7);
stateSpaceColors = [colors1;colors2(2:end,:)];


%%  parameters
lambda = 0.5;
kT = 0.4;
b2 = 0.1;
c1 = 0.1;
c2 = 0.1;
k1 = 1.1;
k2 = 0.95;
ff1 = 0.2; % initial condition of PD-L1
ff2 = 0.1; % initial condition of CCR7
v0 = 1;
yLimits = v0*3;
totalT = 6;
Q = eye(4);

% met params
p1 = 0.5; 
p2 = 0.5;


for treatment_scenario = 1:13 % [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6]

    % determine treatment delays:
    drug1 = [1,1,1,1,1,1];
    drug2 = [1,1,1,1,1,1];
    if (treatment_scenario >= 8)
        drug1(1:(treatment_scenario-7)) = 0;
    end
    
    if (treatment_scenario <=6)
        drug2(1:(7-treatment_scenario)) = 0;
    end
    
        
    
    x0 = [(1 - ff1)*(1 - ff2), ff1*(1-ff2), ff2*(1-ff1), ff1*ff2, v0*(1 - ff1)*(1 - ff2), v0*ff1*(1-ff2), v0*ff2*(1-ff1), v0*ff1*ff2];
    
    % initial met score:
    m0 = p2*( x0(3) + x0(4) )  + p1*( x0(2) + x0(4) );
    met_score = m0;

    %% met score & tumor size
    ttt = [];
    vv = [];
    mm = [];

    % simulate 6 months:
    for i = 1:1:totalT       

        A = payoff(kT, k1*drug1(i), k2*drug2(i), b2, c1, c2);
        K1 = k1*drug1(i);
        K2 = k2*drug2(i);
        
        tplotvec = 0:0.01:1; % 1 month
        [tt, xx]=ode45(@(t,n)rep_ode(t, n, A, Q, lambda), tplotvec, x0,options);
        x0 = xx(end,1:4)./sum(xx(end,1:4));
        x0 = [x0, xx(end,5:8)];

        CCR7 = (xx(:,3) + xx(:,4));
        PDL1 = (xx(:,2) + xx(:,4));

        % calcMetScore returns a
        [delta_met, met_over_time] = calcMetScore(tt,CCR7,PDL1,p1,p2);
        met_over_time = met_score + met_over_time/totalT;
        met_score = met_score + delta_met/totalT;

        %% met score and tumor size
        mm = [mm, met_over_time];
        vv = [vv,sum(xx(:,5:8)')];
        
    end

    figure(1);
    plot(mm,vv,'-','LineWidth',lw,'Color',stateSpaceColors(treatment_scenario,:)); hold on;
    clean_plot(1,'metastatic risk, $\bar{m}$','relative tumor volume, $v$',true);

    % plot start & end (starting dot is always black)
    plot(mm(end),vv(end),'.','MarkerSize',mss,'Color',stateSpaceColors(treatment_scenario,:)); hold on;
    plot(mm(1),vv(1),'.', 'MarkerSize', mss,'Color', 'Black'); hold on;
end


set(gca,'yscale','log');
xlim([0 1]);
ylim([0.05, 10]);























