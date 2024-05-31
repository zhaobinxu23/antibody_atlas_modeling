% parametersetting; 这个模型要解释为什么可能产生慢性感染 chronic infection
clc
clear
%%
mu = -15.5; % 均值
sigma = 0.5; % 标准差
k2 = 1;
k4 = 0.01; %% 0.01
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
        x0(10*(i-1)+j) = prob_A(i)*prob_B(j)*total_B;
    end
end
x0(101) = 1e2;%% virus

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+101) = 0;
    end
end

x0(202) = 1e16; %% 1e16;%% environmental antigen concentration

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


para(1) = 1e-20; 
para(2) = 1e-19;
para(3) = 1e-18; 
para(4) = 1e-17;
para(5) = 1e-16; 
para(6) = 1e-15;
para(7) = 1e-14; 
para(8) = 1e-13;
para(9) = 1e-12; 
para(10) = 1e-11;

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


para(11) = 0.01;%% 0.01
para(12) = 5;
para(13) = 1;
para(14) = 0.5;
para(15) = 1;
para(16) = 0.01*10.5/0.5*1e-16; %% 0.01
para(17) = 10;
para(18) = 1e13;%% 1e13;%% 1e13;



[t y]=ode15s(@pathway_model_many_antibody_immune_res_new,[0 500],x0,[],para,para_new);



for i = 1:100
data_new(i,:) = interp1(t,y(:,i),(0:1:500));
end

data_k_on = zeros(10, 501);

for i = 1:100
    % 计算在 data_k_on 中的行索引
    row_index = fix((i - 1)/10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_on 行
    data_k_on(row_index, :) = data_k_on(row_index, :) + data_new(i, :);
end

data_k_on_index = [-20:1:-11];
average_k_on = sum(data_k_on_index'.*data_k_on)./sum(data_k_on);





data_k_off = zeros(10, 501);

for i = 1:100
    % 计算在data_k_off中的行索引
    row_index = mod(i - 1, 10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_off 行
    data_k_off(row_index, :) = data_k_off(row_index, :) + data_new(i, :);
end

data_k_off_index = [-3:1:6];
average_k_off = sum(data_k_off_index'.*data_k_off)./sum(data_k_off);

data_kd = zeros(19, 501);

for i = 1:100
    % 计算在 data_kd 中的行索引
    row_index = fix((i - 1)/10) - mod(i - 1, 10) + 10;
    
    % 将 data_new 的当前列加到对应的 data_kd 行
    for j = 1:501
        data_kd(row_index, j) = data_kd(row_index, j) + data_new(i, j);
    end
end

data_kd_index = [26:-1:8];
average_kd = sum(data_kd_index'.*data_kd)./sum(data_kd);

data_kd_density = data_kd./sum(data_kd);
%% plot figure 3D
plot(average_k_on);
hold on
plot(average_k_off);
hold on
plot(average_kd);
hold on

%% plot figure 3A
values = [-20:1:-11];
time_points = [0:1:500];
frequencies = data_k_on;
[TimePoints, Values] = meshgrid(time_points, values); % 注意这里的顺序

% 绘制三维频率分布图
figure;
surf(Values, TimePoints, frequencies, 'EdgeColor', 'none');


% 设置 colormap
colormap(parula);

% 添加轴标签和标题
xlabel('Time Points');
ylabel('Values');
zlabel('Frequencies');
title('3D Frequency Distribution');

% 调整视角和视点
view(30, 45);

% 添加网格线
grid on;
%% plot figure 3B
values = [-3:-1:6];
time_points = [0:1:500];
frequencies = data_k_off;
[TimePoints, Values] = meshgrid(time_points, values); % 注意这里的顺序

% 绘制三维频率分布图
figure;
surf(Values, TimePoints, frequencies, 'EdgeColor', 'none');


% 设置 colormap
colormap(parula);

% 添加轴标签和标题
xlabel('Time Points');
ylabel('Values');
zlabel('Frequencies');
title('3D Frequency Distribution');

% 调整视角和视点
view(30, 45);

% 添加网格线
grid on;
%% plot figure 3C
values = [26:-1:8];
time_points = [0:1:500];
frequencies = data_kd;
[TimePoints, Values] = meshgrid(time_points, values); % 注意这里的顺序

% 绘制三维频率分布图
figure;
surf(Values, TimePoints, frequencies, 'EdgeColor', 'none');


% 设置 colormap
colormap(parula);

% 添加轴标签和标题
xlabel('Time Points');
ylabel('Values');
zlabel('Frequencies');
title('3D Frequency Distribution');

% 调整视角和视点
view(30, 45);

% 添加网格线
grid on;

