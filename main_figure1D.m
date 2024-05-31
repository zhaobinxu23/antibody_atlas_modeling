clc
clear
%%
k2 = 1;
k4 = 0.01; %% 0.01
k3 = 0.5;


x0(1) = 1e5;%% Antibody 1
x0(2) = 1e1;%% virus
x0(3) = 0;%% antibody-virus complex
x0(4) = 1e5;%% antibody 2
x0(5) = 0;
x0(6) = 1e5;%% antibody 3
x0(7) = 0;
x0(8) = 1e5;%% antibody 4
x0(9) = 0;

x0(10) = 1e12;%% environmental antigens
x0(11) = x0(1)*k4/(k2-k3);%% E-A complex
x0(12) = x0(4)*k4/(k2-k3);
x0(13) = x0(6)*k4/(k2-k3);
x0(14) = x0(8)*k4/(k2-k3);


para(1) = 1e-12; 
para(2) = 1e-12;
para(3) = 0.9e-12;
para(4) = 0.8e-12;
para(5) = 1e-1;
para(6) = 1;
para(7) = 0.9e-1;
para(8) = 0.8e-1;
para(9) = 0.01;%% 0.01
para(10) = 5;%% feedback constant of virus on antibody regeneration
para(11) = 1;%% virus replication constant
para(12) = 0.5;%% degradtion constant of antibody-antigen complex
para(13) = 1;%% feedback constant of environmental antigens
para(14) = 0.01*10.5/0.5*1e-12; %% k1 for environmental antigens
para(15) = 10;%% k-1 for environmental antigens
para(16) = 4e3;%% environmental replenish constant

% 
% 
% para(2) = 0.01;%% 0.01
% para(3) = 5;%% feedback constant of virus on antibody regeneration
% para(4) = 1;%% virus replication constant
% para(5) = 0.5;%% degradtion constant of antibody-antigen complex
% para(6) = 1;%% feedback constant of environmental antigens
% para(7) = 0.01*10.5/0.5*1e-12; %% k1 for environmental antigens
% para(8) = 10;%% k-1 for environmental antigens
% para(9) = 2e9;%% environmental replenish constant
% para(10) = 1e-1; 



[tt yy]=ode15s(@pathway_model_4_antibody_immune_res_new,[0 500],x0,[],para);

plot(tt,yy(:,1),'linewidth',2);
hold on
plot(tt,yy(:,2),'linewidth',2);
hold on
plot(tt,yy(:,4),'linewidth',2);
hold on
plot(tt,yy(:,6),'linewidth',2);
hold on
plot(tt,yy(:,8),'linewidth',2);
hold on

% dd = interp1(tt,yy(:,1),(0:50:4000));
% plot(dd,'linewidth',2);
% hold on


