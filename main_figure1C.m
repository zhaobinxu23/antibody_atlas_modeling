clc
clear
%%
k2 = 1;
k4 = 0.01; %% 0.01
k3 = 0.5;


x0(1) = 2e11;%% Antibody
x0(2) = 1e1;%% virus
x0(3) = 0;%% antibody-virus complex
x0(4) = 1e12;%% environmental antigens
x0(5) = x0(1)*k4/(k2-k3);%% E-A complex

para(1) = 1e-12; 
para(2) = 0.01;%% 0.01
para(3) = 5;%% feedback constant of virus on antibody regeneration
para(4) = 1;%% virus replication constant
para(5) = 0.5;%% degradtion constant of antibody-antigen complex
para(6) = 1;%% feedback constant of environmental antigens
para(7) = 0.01*10.5/0.5*1e-12; %% k1 for environmental antigens
para(8) = 10;%% k-1 for environmental antigens
para(9) = 2e9;%% environmental replenish constant
para(10) = 1e-1; 



[tt yy]=ode15s(@pathway_model_single_antibody_immune_res,[0 4000],x0,[],para);

plot(tt,yy(:,1),'linewidth',2);
hold on
dd = interp1(tt,yy(:,1),(0:50:4000));
% plot(dd,'linewidth',2);
% hold on


% plot(t,y(:,2),'linewidth',2);
% hold on
% plot(t,y(:,3),'linewidth',2);
% hold on
% plot(t,y(:,4),'linewidth',2);
% hold on
% plot(t,y(:,5),'linewidth',2);
% hold on
%%

k2 = 1;
k4 = 0.01; %% 0.01
k3 = 0.5;


x0(1) = 2e11;%% Antibody
x0(2) = 1e1;%% virus
x0(3) = 0;%% antibody-virus complex
x0(4) = 0;%% environmental antigens
x0(5) = 0;%% E-A complex

para(1) = 1e-12; 
para(2) = 0.01;%% 0.01
para(3) = 5;%% feedback constant of virus on antibody regeneration
para(4) = 1;%% virus replication constant
para(5) = 0.5;%% degradtion constant of antibody-antigen complex
para(6) = 1;%% feedback constant of environmental antigens
para(7) = 0.01*10.5/0.5*1e-12; %% k1 for environmental antigens
para(8) = 10;%% k-1 for environmental antigens
para(9) = 0;%% environmental replenish constant
para(10) = 1e-1; 
para(11) = 0;



[t y]=ode15s(@pathway_model_single_antibody_immune_res_new,[0 4000],x0,[],para);

plot(t,y(:,1),'linewidth',2);
hold on

dd_new = interp1(t,y(:,1),(0:50:4000));

% dd_new = interp1(t,y(:,1),(1000:5:3000));
% plot(dd_new,'linewidth',2);

% plot(t,y(:,2),'linewidth',2);
% hold on
% plot(t,y(:,3),'linewidth',2);
% hold on
% plot(t,y(:,4),'linewidth',2);
% hold on
% plot(t,y(:,5),'linewidth',2);