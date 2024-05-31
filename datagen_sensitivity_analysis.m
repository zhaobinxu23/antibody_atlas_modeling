function [integral_value,max_virus,max_symptom,max_elisa,final_elisa] = datagen_sensitivity_analysis(theta)

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


mu_2 = theta(3); % 均值
sigma_2 = theta(4); % 标准差


prob_B(1) = normcdf(theta(3)-4, mu_2, sigma_2);
% 计算累积分布概率
for i = 2:5
prob_B(i) = normcdf(theta(3)-4+i-1, mu_2, sigma_2) - normcdf(theta(3)-4+i-2, mu_2, sigma_2);
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

para_new(1) = 10^(theta(3)-4.5); 
para_new(2) = 10^(theta(3)-3.5); 
para_new(3) = 10^(theta(3)-2.5); 
para_new(4) = 10^(theta(3)-1.5); 
para_new(5) = 10^(theta(3)-0.5);  
para_new(6) = 10^(theta(3)+0.5); 
para_new(7) = 10^(theta(3)+1.5);  
para_new(8) = 10^(theta(3)+2.5); 
para_new(9) = 10^(theta(3)+3.5);  
para_new(10) = 10^(theta(3)+4.5); 


para(11) = 0.01;
para(12) = theta(5);%% 5
para(13) = theta(6);%% 1
para(14) = theta(7);%% 0.5
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

[ttt,yyy]=ode15s(@pathway_model_many_antibody_immune_res_new,[0 200],x0,[],para,para_new);
[t_final,y_final]=ode15s(@pathway_model_many_antibody_immune_res_new,[0 10000],x0,[],para,para_new);

integral_value = trapz(ttt, yyy(:,101));
max_virus = max(yyy(:,101));

for i = 1:100
data_virus_antibody_complex(i,:) = interp1(ttt,yyy(:,i+101),(0:1:200));
end
data_virus_antibody_complex_overall = sum(data_virus_antibody_complex);

max_symptom = max(data_virus_antibody_complex_overall);



for i = 1:100
data_new_2(i,:) = interp1(ttt,yyy(:,i),(0:1:200));
end

timepoints = (0:1:200);

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
        x0(10*(i-1)+j+202) = 0;
    end
end

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

para_new(1) = 10^(theta(3)-4.5); 
para_new(2) = 10^(theta(3)-3.5);
para_new(3) = 10^(theta(3)-2.5); 
para_new(4) = 10^(theta(3)-1.5);
para_new(5) = 10^(theta(3)-0.5); 
para_new(6) = 10^(theta(3)+0.5);
para_new(7) = 10^(theta(3)+1.5); 
para_new(8) = 10^(theta(3)+2.5);
para_new(9) = 10^(theta(3)+3.5); 
para_new(10) = 10^(theta(3)+4.5);


para_2(11) = 0;%% 0.01
para_2(12) = 0;%% 5
para_2(13) = 0;%% 1
para_2(14) = 0;%% 0.5
para_2(15) = 1;
para_2(16) = 0.01*10.5/0.5*1e-16;
para_2(17) = 10;
para_2(18) = 0;



[t, zz]=ode15s(@pathway_model_many_antibody_immune_res_new,[0 10],x0,[],para_2,para_new);

out_elisa(k) = 1e16-interp1(t,zz(:,101),10);
end

max_elisa = max(out_elisa);



for i = 1:100
    x0(i) =  interp1(t_final,y_final(:,i),(10000));
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
        x0(10*(i-1)+j+202) = 0;
    end
end



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

para_new(1) = 10^(theta(3)-4.5); 
para_new(2) = 10^(theta(3)-3.5);
para_new(3) = 10^(theta(3)-2.5); 
para_new(4) = 10^(theta(3)-1.5);
para_new(5) = 10^(theta(3)-0.5); 
para_new(6) = 10^(theta(3)+0.5);
para_new(7) = 10^(theta(3)+1.5); 
para_new(8) = 10^(theta(3)+2.5);
para_new(9) = 10^(theta(3)+3.5); 
para_new(10) = 10^(theta(3)+4.5);


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

[t_new_new, zz]=ode15s(@pathway_model_many_antibody_immune_res_new,[0 10],x0,[],para_2,para_new);

final_elisa = 1e16-interp1(t_new_new,zz(:,101),10);
end



