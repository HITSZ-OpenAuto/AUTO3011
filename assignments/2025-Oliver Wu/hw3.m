clc;clear;close all;

% 读取整个文件到矩阵
prefix = "001";
data = load("motion control/lab3/"+prefix+".fwtxt");
t_x = (0:length(data(:,1))-1).*0.002;

% 将每一列分配给不同的数组
col1 = data(:, 1);  % 第一列，规划x
col2 = data(:, 2);  % 第二列，规划y
col3 = data(:, 3);  % 第三列，编码器x
col4 = data(:, 4);  % 第四列，编码器y

% 初始化变量
num_points = length(col1);
contour_errors = zeros(num_points, 1);  % 存储轮廓误差
nearest_points = zeros(num_points, 2);  % 存储最近点的坐标 [x, y]
segment_ids = zeros(num_points, 1);  % 存储各点对应最近点所在段

% 第一个点的处理
current_segment_start = 1;  % 当前搜索线段的起始索引
nearest_points(1, :) = [col1(1), col2(1)];  % 第一个最近点就是第一个规划点
contour_errors(1) = sqrt((col3(1) - col1(1))^2 + (col4(1) - col2(1))^2);
segment_ids(1) = 1;

% 遍历所有点（从第二个点开始）
for i = 2:num_points
    current_point = [col3(i), col4(i)];  % 当前实际位置
    
    % 搜索范围：从上一个最近点所在的线段开始，到当前规划点
    min_dist = inf;
    best_nearest_point = [0, 0];
    best_segment_start = current_segment_start;
    
    % 在可能的线段中寻找最近点
    for seg_start = current_segment_start:min(i-1, num_points-1)
        seg_end = seg_start + 1;
        if seg_start==current_segment_start
            % 第一条线段要从上个最近点开始算
            segment_start_point = nearest_points(i-1, :);
        else
            segment_start_point = [col1(seg_start), col2(seg_start)];
        end
        segment_end_point = [col1(seg_end), col2(seg_end)];
        
        [dist, nearest_pt] = point_to_segment_distance(current_point, ...
                                                      segment_start_point, ...
                                                      segment_end_point);
        
        if dist < min_dist
            min_dist = dist;
            best_nearest_point = nearest_pt;
            best_segment_start = seg_start;
        end
    end
    
    % 更新结果
    contour_errors(i) = min_dist;
    nearest_points(i, :) = best_nearest_point;
    segment_ids(i) = best_segment_start;
    current_segment_start = best_segment_start;
end

% 绘制结果
figure();
% 图1：轨迹和轮廓误差可视化
plot(col1, col2, 'b-', 'LineWidth', 2, 'DisplayName', '规划轨迹');
hold on;
plot(col3, col4, 'r-', 'LineWidth', 1, 'DisplayName', '实际轨迹');
% plot(nearest_points(:,1), nearest_points(:,2), 'g--', 'LineWidth', 1, 'DisplayName', '最近点轨迹');
legend();

% 绘制一些误差连线示例（每隔若干个点画一条）
for i = 1:100:num_points
    h = plot([col3(i), nearest_points(i,1)], [col4(i), nearest_points(i,2)], ...
         'm-', 'LineWidth', 0.5, 'DisplayName', '实际点与轮廓上最近点连线');
    if i ~= 1
        h.HandleVisibility = 'off';  % 不出现在图例中
    end
end

xlabel('X坐标');
ylabel('Y坐标');
title('轨迹和轮廓误差');
grid on;
axis equal;

figure;

% 图2.1：轮廓误差随时间变化
subplot(2,1,1);
plot(t_x, contour_errors, 'k-', 'LineWidth', 1.5); % 转换为毫米
xlabel('时间 (s)');
ylabel('轮廓误差 (mm)');
title('轮廓误差时间序列');
grid on;

% 图2.2：轮廓误差统计
subplot(2,1,2);
histogram(contour_errors, 50, 'FaceColor', 'cyan', 'EdgeColor', 'blue');
xlabel('轮廓误差 (mm)');
ylabel('频数');
title('轮廓误差分布直方图');
grid on;

% 保存结果
results = table(t_x', col3, col4, nearest_points(:,1), nearest_points(:,2), ...
               contour_errors, 'VariableNames', ...
               {'Time', 'Actual_X', 'Actual_Y', 'Nearest_X', 'Nearest_Y', 'Contour_Error_mm'});
% save("contour_error_"+prefix+".mat", "results", '-mat');
writetable(results, "contour_error_results_"+prefix+".csv");

% 显示主要统计结果
fprintf('\n=== 轮廓误差分析结果 ===\n');
fprintf('最大轮廓误差: %.3f mm\n', max(contour_errors));
fprintf('平均轮廓误差: %.3f mm\n', mean(contour_errors));
fprintf('轮廓误差均方根: %.3f mm\n', rms(contour_errors));

% 计算点到线段的距离函数
function [min_distance, nearest_point] = point_to_segment_distance(point, seg_start, seg_end)
    % 计算点到线段的最短距离和最近点
    v = seg_end - seg_start;
    w = point - seg_start;
    
    c1 = dot(w, v);
    if c1 <= 0
        % 最近点是线段起点
        nearest_point = seg_start;
        min_distance = norm(point - seg_start);
        return;
    end
    
    c2 = dot(v, v);
    if c2 <= c1
        % 最近点是线段终点
        nearest_point = seg_end;
        min_distance = norm(point - seg_end);
        return;
    end
    
    % 最近点在线段内部
    b = c1 / c2;
    nearest_point = seg_start + b * v;
    min_distance = norm(point - nearest_point);
end