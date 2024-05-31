% parametersetting; 这个模型要解释为什么可能产生慢性感染 chronic infection
%% run Main_immune_imprinting_many_antibodies first

% for i = 1:100
% data_new(i,:) = interp1(t,y(:,i),(0:10:1000));
% end
preinfection_1;
%% ELISA results


mu = -15.5; % 均值 -15.5
sigma = 0.5; % 标准差
k2 = 1;
k4 = 0.01;
k3 = 0.5;
k_1 = 1;


prob_A(1) = normcdf(-19.5, mu, sigma);
% 计算累积分布概率
for i = 2:5
prob_A(i) = normcdf(-19.5+i-1, mu, sigma) - normcdf(-19.5+i-2, mu, sigma);
end
for i = 6:10
prob_A(i) = prob_A(11-i);
end


mu = 1.5; % 均值
sigma = 0.8; % 标准差


prob_B(1) = normcdf(-2.5, mu, sigma);
% 计算累积分布概率
for i = 2:5
prob_B(i) = normcdf(-2.5+i-1, mu, sigma) - normcdf(-2.5+i-2, mu, sigma);
end
for i = 6:10
prob_B(i) = prob_B(11-i);
end

total_B = 1e15;

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j) = prob_A(i)*prob_B(j)*total_B*(1-1e-10)+1e-10*data_new(10*(i-1)+j,1001);
        zz0(10*(i-1)+j) = prob_A(i)*prob_B(j)*total_B*(1-1e-10);
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



para(1) = 10^(-20); 
para(2) = 10^(-19);
para(3) = 10^(-18); 
para(4) = 10^(-17);
para(5) = 10^(-16); 
para(6) = 10^(-15);
para(7) = 10^(-14); 
para(8) = 10^(-13);
para(9) = 10^(-12); 
para(10) = 10^(-11);

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
para(12) = 5;
para(13) = 1;
para(14) = 0.5;
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

[t_new z]=ode15s(@pathway_model_many_antibody_immune_res_new,[0 1000],x0,[],para,para_new);



for i = 1:100
data_new_sin(i,:) = interp1(t_new,z(:,i),(0:1:1000));
end
% 
% plot((200:1:400),sum(data_new_sin),'linewidth',2);
% hold on


data_k_on_sin = zeros(10, 1001);

for i = 1:100
    % 计算在 data_k_on 中的行索引
    row_index = fix((i - 1)/10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_on 行
    data_k_on_sin(row_index, :) = data_k_on_sin(row_index, :) + data_new_sin(i, :);
end

data_k_off_sin = zeros(10, 1001);

for i = 1:100
    % 计算在data_k_off中的行索引
    row_index = mod(i - 1, 10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_off 行
    data_k_off_sin(row_index, :) = data_k_off_sin(row_index, :) + data_new_sin(i, :);
end

data_kd_sin = zeros(19, 1001);

for i = 1:100
    % 计算在 data_kd 中的行索引
    row_index = fix((i - 1)/10) - mod(i - 1, 10) + 10;
    
    % 将 data_new 的当前列加到对应的 data_kd 行
    for j = 1:1001
        data_kd_sin(row_index, j) = data_kd_sin(row_index, j) + data_new_sin(i, j);
    end
end

% values = [26:-1:8];
% time_points = [0:1:1000];
% frequencies = data_kd_sin;
% [TimePoints, Values] = meshgrid(time_points, values); % 注意这里的顺序
% 
% % 绘制三维频率分布图
% figure;
% surf(Values, TimePoints, frequencies, 'EdgeColor', 'none');
% 
% 
% % 设置 colormap
% colormap(parula);
% 
% % 添加轴标签和标题
% xlabel('Time Points');
% ylabel('Values');
% zlabel('Frequencies');
% title('3D Frequency Distribution');
% 
% 
% % 调整视角和视点
% view(30, 45);
% 
% % 添加网格线
% grid on;

