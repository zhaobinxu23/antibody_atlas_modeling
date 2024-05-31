
antigenic_sin_atlas_updated_version_1
for i = 1:100
    plot(t, y(:,i), 'linewidth', 2, 'color', 'b');
    hold on
    plot(t_new, z(:,i), 'linewidth', 2, 'color', 'r');
    
    xlabel('Time', 'FontWeight', 'bold');
    ylabel('Concentration', 'FontWeight', 'bold');
    title(['Plot-antibody', num2str(i), '-dynamics']);
    
    legend('without antigenic sin', 'with antigenic sin', 'Location', 'best');

    % 可以在此处添加其他修饰图形的代码，例如添加网格、调整坐标轴范围等
    % ...

    % 保存每个图像为不同的文件（可选）
    filename = ['plot', num2str(i), '.png'];
    saveas(gcf, filename);
    close all
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