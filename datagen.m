function [out] = datagen(timepoints, theta)
mu = theta(1); % 均值-15.5
sigma = theta(2); % 标准差 0.5
k2 = 1;
k4 = 0.01;
k3 = 0.5;
k_1 = 1;


prob_A(1) = normcdf(theta(1)-4, mu, sigma);
% 计算累积分布概率
for i = 2:5
prob_A(i) = normcdf(theta(1)-4+i-1, mu, sigma) - normcdf(theta(1)-4+i-2, mu, sigma);
end
for i = 6:10
prob_A(i) = prob_A(11-i);
end


mu_2 = 1.5; % 均值
sigma_2 = 0.8; % 标准差


prob_B(1) = normcdf(-2.5, mu_2, sigma_2);
% 计算累积分布概率
for i = 2:5
prob_B(i) = normcdf(-2.5+i-1, mu_2, sigma_2) - normcdf(-2.5+i-2, mu_2, sigma_2);
end
for i = 6:10
prob_B(i) = prob_B(11-i);
end

total_B = 1e15;

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j) = prob_A(i)*prob_B(j)*total_B;
    end
end
x0(101) = 1e2;%% virus

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+101) = 0;
    end
end

x0(202) = 1e16;%% environmental antigen concentration

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+202) = x0(10*(i-1)+j)*k4/(k2-k3);
    end
end

% x0(1) = 1e8;% antibody 1
% x0(2) = 1e14;% antibody 2
% x0(3) = 1e2;% virus
% x0(4) = 0;% antibody-1-virus complex
% x0(5) = 0;% antibody-2-virus complex


para(1) = 10^(theta(1)-4.5); 
para(2) = 10^(theta(1)-3.5); 
para(3) = 10^(theta(1)-2.5); 
para(4) = 10^(theta(1)-1.5); 
para(5) = 10^(theta(1)-0.5);  
para(6) = 10^(theta(1)+0.5); 
para(7) = 10^(theta(1)+1.5);  
para(8) = 10^(theta(1)+2.5); 
para(9) = 10^(theta(1)+3.5);  
para(10) = 10^(theta(1)+4.5); 

para_new(1) = 1e-3; 
para_new(2) = 1e-2;
para_new(3) = 1e-1; 
para_new(4) = 1e0;
para_new(5) = 1e1; 
para_new(6) = 1e2;
para_new(7) = 1e3; 
para_new(8) = 1e4;
para_new(9) = 1e5; 
para_new(10) = 1e6;


para(11) = 0.01;
para(12) = theta(3);%% 5
para(13) = theta(4);%% 1
para(14) = theta(5);%% 0.5
para(15) = 1;
para(16) = 0.01*10.5/0.5*1e-16;
para(17) = 10;
para(18) = 1e13;

% para(2) = 1e-15; 
% para(3) = 1; 
% para(4) = 0.01; 
% para(5) = 5;
% para(6) = 1;
% para(7) = 0.5;
% para(8) = 1;

[t,y]=ode15s(@pathway_model_many_antibody_immune_res_new,[0 10000],x0,[],para,para_new);

for i = 1:100
data_new_2(i,:) = interp1(t,y(:,i),(0:10:10000));
end

for k = 1:length(timepoints)


for i = 1:100
    x0(i) = data_new_2(i,timepoints(k)+1);
end

x0(101) = 1e16;%% virus

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+101) = 0;
    end
end

x0(202) = 0;%% environmental antigen concentration

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+202) = x0(10*(i-1)+j)*k4/(k2-k3);
    end
end

% x0(1) = 1e8;% antibody 1
% x0(2) = 1e14;% antibody 2
% x0(3) = 1e2;% virus
% x0(4) = 0;% antibody-1-virus complex
% x0(5) = 0;% antibody-2-virus complex


para_2(1) = 10^(theta(1)-4.5); 
para_2(2) = 10^(theta(1)-3.5); 
para_2(3) = 10^(theta(1)-2.5); 
para_2(4) = 10^(theta(1)-1.5); 
para_2(5) = 10^(theta(1)-0.5);  
para_2(6) = 10^(theta(1)+0.5); 
para_2(7) = 10^(theta(1)+1.5);  
para_2(8) = 10^(theta(1)+2.5); 
para_2(9) = 10^(theta(1)+3.5);  
para_2(10) = 10^(theta(1)+4.5); 

para_new(1) = 1e-3; 
para_new(2) = 1e-2;
para_new(3) = 1e-1; 
para_new(4) = 1e0;
para_new(5) = 1e1; 
para_new(6) = 1e2;
para_new(7) = 1e3; 
para_new(8) = 1e4;
para_new(9) = 1e5; 
para_new(10) = 1e6;


para_2(11) = 0;%% 0.01
para_2(12) = 0;%% 5
para_2(13) = 0;%% 1
para_2(14) = 0;%% 0.5
para_2(15) = 1;
para_2(16) = 0.01*10.5/0.5*1e-16;
para_2(17) = 10;
para_2(18) = 0;

% para(2) = 1e-15; 
% para(3) = 1; 
% para(4) = 0.01; 
% para(5) = 5;
% para(6) = 1;
% para(7) = 0.5;
% para(8) = 1;

[t, zz]=ode15s(@pathway_model_many_antibody_immune_res_new,[0 10],x0,[],para_2,para_new);

out(k) = 1e16-interp1(t,zz(:,101),10);

end



end