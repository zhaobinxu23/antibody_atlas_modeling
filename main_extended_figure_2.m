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



[t y]=ode15s(@pathway_model_many_antibody_immune_res_new,[0 1000],x0,[],para,para_new);


for i = 1:100
data_new(i,:) = interp1(t,y(:,i),(0:1:1000));
end

data_k_on = zeros(10, 1001);

for i = 1:100
    % 计算在 data_k_on 中的行索引
    row_index = fix((i - 1)/10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_on 行
    data_k_on(row_index, :) = data_k_on(row_index, :) + data_new(i, :);
end




data_k_off = zeros(10, 1001);

for i = 1:100
    % 计算在data_k_off中的行索引
    row_index = mod(i - 1, 10) + 1;
    
    % 将 data_new 的当前列加到对应的 data_k_off 行
    data_k_off(row_index, :) = data_k_off(row_index, :) + data_new(i, :);
end

data_kd = zeros(19, 1001);

for i = 1:100
    % 计算在 data_kd 中的行索引
    row_index = fix((i - 1)/10) - mod(i - 1, 10) + 10;
    
    % 将 data_new 的当前列加到对应的 data_kd 行
    for j = 1:1001
        data_kd(row_index, j) = data_kd(row_index, j) + data_new(i, j);
    end
end
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

x0(202) = 0;%% environmental antigen concentration

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+202) = 0;
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
para(18) = 0;%% 1e13;



[t y]=ode15s(@pathway_model_many_antibody_immune_res_new,[0 500],x0,[],para,para_new);

for i = 1:100
    plot(t, y(:,i), 'linewidth', 2);
    xlabel('Time', 'FontWeight', 'bold');
    ylabel('Concentration', 'FontWeight', 'bold');
    title(['Plot-antibody- ', num2str(i),'-dynamics']);
    % 可以在此处添加其他修饰图形的代码，例如添加网格、调整坐标轴范围等
    % ...
    % 保存每个图像为不同的文件（可选）
    filename = ['plot', num2str(i), '.png'];
    saveas(gcf, filename);
end



fig = figure('Position', [0 0 3000 2000], 'PaperUnits', 'inches', 'PaperPosition', [0 0 30 20], 'PaperSize', [30 20]);

% 分成10行10列，并设置间距
for i = 1:100
subplot(10,10,i);
% 读取每个图像文件
filename = ['plot', num2str(i), '.png'];
img = imread(filename);
% 在subplot中显示图像
imagesc(img);
axis off;
end

saveas(fig, 'merged_image.png');