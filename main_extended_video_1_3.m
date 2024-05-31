% parametersetting; 这个模型要解释为什么可能产生慢性感染 chronic infection
clc
clear
%%
mu = -15.5; % 均值
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


para(11) = 0.01;
para(12) = 5;
para(13) = 1;
para(14) = 0.5;
para(15) = 1;
para(16) = 0.01*10.5/0.5*1e-16;
para(17) = 10;
para(18) = 1e13;



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




data_k_off = zeros(10, 501);

for i = 1:100
    % 计算在data_k_off中的行索引
    row_index = mod(i - 1, 10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_off 行
    data_k_off(row_index, :) = data_k_off(row_index, :) + data_new(i, :);
end

data_kd = zeros(19, 501);

for i = 1:100
    % 计算在 data_kd 中的行索引
    row_index = fix((i - 1)/10) - mod(i - 1, 10) + 10;
    
    % 将 data_new 的当前列加到对应的 data_kd 行
    for j = 1:501
        data_kd(row_index, j) = data_kd(row_index, j) + data_new(i, j);
    end
end




timePoints = [0:1:500];

% 假设你有一个名为 data 的绘图数据，其中每一列代表一个时间点的数据
% 假设 data 是一个 10x5 的矩阵，表示有 10 个数据点，5 个时间点
data = data_kd;

% 假设你有另一个名为 data2 的绘图数据，与 data 具有相同的维度
% data2 = data_kd_no_self;

% 初始化 Y 轴范围的最大和最小值
yMin = min([data(:)]);
yMax = max([data(:)]);

% 创建一个新的图形窗口
figure;

% 创建视频写入对象
writerObj = VideoWriter('extended_video_3.avi');
open(writerObj);

% 定义X轴坐标映射关系
x_map = [26:-1:8]; % 1对应26，2对应25，...，20对应8
% x_map = [0:1:9]; % 1对应26，2对应25，...，20对应8
% 绘制每个时间点的图像，并将帧写入视频文件
for tt = 1:length(timePoints)
    % 绘制曲线 data
    plot(x_map, data(:, tt), 'LineWidth', 2, 'Color', [0.2, 0.6, 0.8]); % 使用深蓝色绘制曲线

    
    % 计算阴影区域
%     x = x_map;
%     y1 = data(:, tt);
%  
%     fill([x, fliplr(x)], [y1',fliplr(y1')], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % 设置 Y 轴范围
    ylim([yMin, yMax]);
    
    hold off;
    
    title(['Time Point: ', num2str(timePoints(tt))], 'FontSize', 14); % 更大的标题字体
    xlabel('ln(Kd)', 'FontSize', 12); % X 轴标签
    ylabel('Concentration', 'FontSize', 12); % Y 轴标签
    
    % 捕获当前图像帧
    frame = getframe(gcf);
    
    % 将当前帧写入视频文件
    writeVideo(writerObj, frame);
end

% 关闭视频写入对象
close(writerObj);
