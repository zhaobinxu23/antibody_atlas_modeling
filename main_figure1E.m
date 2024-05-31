% parametersetting; 这个模型要解释为什么可能产生慢性感染 chronic infection
clc
clear

x0(1) = 1e14;% antibody
x0(2) = 1e2;% virus
x0(3) = 0;% complex

para(1) = 1e-15;
para(2) = 1e1;
para(3) = 5;
para(4) = 1;
para(5) = 0.5;


[t y]=ode15s(@pathway_model_one_antibody_immune_res,[0 100],x0,[],para);

% plot(t,y(:,1),'linewidth',2);
% hold on
plot(t,y(:,2),'linewidth',2);
hold on
plot(t,y(:,3),'linewidth',2);
hold on
%% 
x0(1) = 1e14;
x0(2) = 1e2;
x0(3) = 0;

para(1) = 2e-15;
para(2) = 1e1;
para(3) = 5;
para(4) = 1;
para(5) = 0.5;


[t y]=ode15s(@pathway_model_one_antibody_immune_res,[0 100],x0,[],para);

% plot(t,y(:,1),'linewidth',2);
% hold on
plot(t,y(:,2),'linewidth',2);
hold on
plot(t,y(:,3),'linewidth',2);
hold on



