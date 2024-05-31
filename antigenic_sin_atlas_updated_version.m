% parametersetting; 这个模型要解释为什么可能产生慢性感染 chronic infection
%% run Main_immune_imprinting_many_antibodies first

% for i = 1:100
% data_new(i,:) = interp1(t,y(:,i),(0:10:1000));
% end
pre_infection_antigenic_sin;
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
        x0(10*(i-1)+j) = prob_A(i)*prob_B(j)*total_B*(1-1e-14)+1e-14*data_new(10*(i-1)+j,201);
        zz0(10*(i-1)+j) = prob_A(i)*prob_B(j)*total_B*(1-1e-14);
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

[t_new z]=ode15s(@pathway_model_many_antibody_immune_res_new,[200 400],x0,[],para,para_new);

for i = 1:200

strong_antibody_new(i) = interp1(t_new,z(:,51),200+i)+interp1(t_new,z(:,61),200+i)+interp1(t_new,z(:,62),200+i)+interp1(t_new,z(:,71),200+i)+interp1(t_new,z(:,72),200+i)...
    +interp1(t_new,z(:,73),200+i)+interp1(t_new,z(:,81),200+i)+interp1(t_new,z(:,82),200+i)+interp1(t_new,z(:,83),200+i)+interp1(t_new,z(:,84),200+i)+interp1(t_new,z(:,91),200+i)...
    +interp1(t_new,z(:,92),200+i)+interp1(t_new,z(:,93),200+i)+interp1(t_new,z(:,94),200+i)+interp1(t_new,z(:,95),200+i);
end

for i = 1:200
    cross_over_1(i) = (interp1(t_new,z(:,51),200+i)*(x0(51)-zz0(51))/x0(51) + interp1(t_new,z(:,61),200+i)*(x0(61)-zz0(61))/x0(61)+interp1(t_new,z(:,62),200+i)*(x0(62)-zz0(62))/x0(62)...
        +interp1(t_new,z(:,71),200+i)*(x0(71)-zz0(71))/x0(71)+interp1(t_new,z(:,72),200+i)*(x0(72)-zz0(72))/x0(72)+interp1(t_new,z(:,73),200+i)*(x0(73)-zz0(73))/x0(73)...
        +interp1(t_new,z(:,81),200+i)*(x0(81)-zz0(81))/x0(81)+interp1(t_new,z(:,82),200+i)*(x0(82)-zz0(82))/x0(82)+interp1(t_new,z(:,83),200+i)*(x0(83)-zz0(83))/x0(83)...
        +interp1(t_new,z(:,84),200+i)*(x0(84)-zz0(84))/x0(84)+interp1(t_new,z(:,91),200+i)*(x0(91)-zz0(91))/x0(91)+interp1(t_new,z(:,92),200+i)*(x0(92)-zz0(92))/x0(92)...
        +interp1(t_new,z(:,93),200+i)*(x0(93)-zz0(93))/x0(93)+interp1(t_new,z(:,94),200+i)*(x0(94)-zz0(94))/x0(94)+interp1(t_new,z(:,95),200+i)*(x0(95)-zz0(95))/x0(95))/strong_antibody_new(i);
    cross_over_2(i) = cross_over_1(i)*strong_antibody_new(i)/strong_antibody;
end



% plot(t_new,z(:,101),'linewidth',2);
% hold on
for i = 1:100
data_new_sin(i,:) = interp1(t_new,z(:,i+101),(200:1:400));
end
% 
plot((200:1:400),sum(data_new_sin),'linewidth',2);
hold on


data_k_on_sin = zeros(10, 201);

for i = 1:100
    % 计算在 data_k_on 中的行索引
    row_index = fix((i - 1)/10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_on 行
    data_k_on_sin(row_index, :) = data_k_on_sin(row_index, :) + data_new_sin(i, :);
end

data_k_off_sin = zeros(10, 201);

for i = 1:100
    % 计算在data_k_off中的行索引
    row_index = mod(i - 1, 10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_off 行
    data_k_off_sin(row_index, :) = data_k_off_sin(row_index, :) + data_new_sin(i, :);
end

data_kd_sin = zeros(19, 201);

for i = 1:100
    % 计算在 data_kd 中的行索引
    row_index = fix((i - 1)/10) - mod(i - 1, 10) + 10;
    
    % 将 data_new 的当前列加到对应的 data_kd 行
    for j = 1:201
        data_kd_sin(row_index, j) = data_kd_sin(row_index, j) + data_new_sin(i, j);
    end
end

% values = [-26:1:-8];
% time_points = [0:5:1000];
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

