clc;clear;close all;

% 读取整个文件到矩阵
data = load('EppToothcurve.txt'); 

% 将每一列分配给不同的数组
col1 = data(:, 1);
col2 = data(:, 2);
% 第一步：绘制出这个图形
axis equal;
plot(col1, col2);
% scatter(col1, col2);

% 第二步：分区域
seg = segment_curve(col1, col2);

% 第三步：绘制图形
plot_segmentation_results(col1, col2, seg);

function segments = segment_curve(x, y)
    % 主要分段算法
    % 输入：x, y - 曲线坐标
    % 输出：segments - 分段信息
    
    n = length(x);
    segments = [];
    current_segment_start = 1;
    
    % 参数设置
    angle_threshold = 45; % 尖锐转角阈值(度)
    line_length_threshold = 1.0; % 长线段长度阈值
    line_fit_threshold = 0.02; % 直线拟合误差阈值
    
    i = 1;
    while i < n - 2
        % 检测尖锐转角
        if detect_sharp_corner(x, y, i, angle_threshold)
            % 找到分段点
            segment_end = i;
            segments = add_segment(segments, i, i+2, 3);
            
            % 添加分段
            if segment_end >= current_segment_start + 2  % 至少三个点
                segments = detect_long_segments(x, y, ...
                    current_segment_start, segment_end, ...
                    line_length_threshold, line_fit_threshold, ...
                    segments);
                current_segment_start = segment_end+2;
                i = current_segment_start;
            else
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end
    
    % 添加最后一段
    if current_segment_start < n
        segments = detect_long_segments(x, y, ...
                    current_segment_start, n, ...
                    line_length_threshold, line_fit_threshold, ...
                    segments);
    end

    segments = merge_short_segments(segments, 5, 15);
end

function is_sharp = detect_sharp_corner(x, y, idx, threshold_deg)
    % 检测尖锐转角, 使用相邻三点计算转角
    
    if idx > length(x)-2
        is_sharp = false;
        return;
    end
    
    % 计算向量
    v1 = [x(idx+1) - x(idx), y(idx+1) - y(idx)];
    v2 = [x(idx+2) - x(idx+1), y(idx+2) - y(idx+1)];
    
    % 计算夹角
    cos_angle = dot(v1, v2) / (norm(v1) * norm(v2));
    cos_angle = max(min(cos_angle, 1), -1); % 确保在有效范围内
    angle_deg = acosd(cos_angle);

    is_sharp = (min(angle_deg, 360-angle_deg) > threshold_deg);
end

function segments = add_segment(segments, start_idx, end_idx, type, linear_params)
    % 添加新分段
    if type == 1  % 直线
        new_segment = struct('start', start_idx, 'end', end_idx, 'type', 1, 'linear_params', linear_params);
    elseif type == 2  % 曲线
        new_segment = struct('start', start_idx, 'end', end_idx, 'type', 2, 'linear_params', []);
    elseif type == 3  % 急转
        new_segment = struct('start', start_idx, 'end', end_idx, 'type', 3, 'linear_params', []);
    else  % 错误类型
        error("错误的分段类型！");
    end
    segments = [segments, new_segment];
end

function segments = detect_long_segments(x, y, start_idx, end_idx, length_threshold, fit_threshold, segments)
    % 检测长线段区域
    % 点数太少
    if (end_idx - start_idx) <= 2
        [linear_fit, fit_error] = check_linear_fit(x(start_idx:end_idx), ...
                                                   y(start_idx:end_idx));
        if fit_error < fit_threshold
            segments = add_segment(segments, start_idx, end_idx, 1, linear_fit);
        else
            segments = add_segment(segments, start_idx, end_idx, 2);
        end
        return;
    end
    % 直线拟合检查
    current_start = start_idx;
    current_end = start_idx + 2;
    last_linear_fit = [];
    last_fit_error = Inf;
    while current_end <= end_idx
        [linear_fit, fit_error] = check_linear_fit(x(current_start:current_end), ...
                                                  y(current_start:current_end));
        % 若拟合误差很小, 或刚开始计算
        if ((fit_error < fit_threshold)) || (current_end-current_start <= 2)
            last_linear_fit = linear_fit;
            last_fit_error = fit_error;
            if end_idx == current_end
                segments = add_segment(segments, current_start, current_end, 1, linear_fit);
            end
            current_end = current_end+1;
        else
            segment_length = sqrt((x(current_end) - x(current_start))^2 + ...
                                  (y(current_end) - y(current_start))^2);
            % 较长（长直线）判断
            if segment_length > length_threshold
                segments = add_segment(segments, current_start, current_end-1, 1, last_linear_fit);
            else
                segments = add_segment(segments, current_start, current_end-1, 2);
            end
            last_fit_error = Inf;
            last_linear_fit = [];
            if end_idx - (current_end-1) <= 2  % 最后阶段
                segments = add_segment(segments, current_end-1, end_idx, 2);
                break;
            else
                current_start = current_end-1;
                current_end = current_end + 1;
            end
        end
    end
end

% 连直线算弦差
function [linear_params, fit_error] = check_linear_fit(x_seg, y_seg)
    % 检查段是否可以很好地用直线拟合
    
    n = length(x_seg);
    if n < 3
        linear_params = [];
        fit_error = inf;
        return;
    end
    [slope, inter] = line_parameters(x_seg(1), y_seg(1), x_seg(end), y_seg(end));
    linear_params = [slope; inter];
    % 计算误差
    fit_error = 0;
    A = [x_seg(1), y_seg(1)];
    B = [x_seg(end), y_seg(end)];
    for i=2:n-1
        P = [x_seg(i), y_seg(i)];
        error = point_to_seg_dist(P,A,B);
        fit_error = max(error, fit_error);
    end
end

function d = point_to_seg_dist(P, A, B)
    AB = B - A;  ab2 = dot(AB,AB);
    if ab2 < eps, d = norm(P-A); return; end
    t = dot(P - A, AB) / ab2;
    t = max(0,min(1,t));
    Q = A + t*AB;
    d = norm(P - Q);
end

function [slope, intercept] = line_parameters(x1, y1, x2, y2)
% 计算通过两点的直线参数
% 输入: 两点坐标 (x1,y1), (x2,y2)
% 输出: 斜率slope, 截距intercept

    if x1 == x2
        slope=Inf;
        intercept=Inf;
        return;
    end
    
    slope = (y2 - y1) / (x2 - x1);
    intercept = y1 - slope * x1;
    
    % 可选：显示结果
%     fprintf('直线方程: y = %.4fx + %.4f\n', slope, intercept);
end

function segments = merge_short_segments(segments, min_points, max_points)
    % 合并过短的段
    
    if length(segments) <= 1
        return;
    end
    
    i = 2;
    while i <= length(segments)
        % 只合并曲线
        if segments(i-1).type == 1 || segments(i).type == 1
            i = i+1;
            continue;
        end
        prev_length = segments(i-1).end - segments(i-1).start + 1;
        curr_length = segments(i).end - segments(i).start + 1;
        % 合并的段不应太长，以利于插值
        if prev_length > max_points
            i = i+1;
            continue;
        end
        
        if prev_length < min_points || curr_length < min_points
            % 合并段
            segments(i-1).end = segments(i).end;
            segments(i) = [];
        else
            i = i + 1;
        end
    end
end

function plot_segmentation_results(x, y, segments)
    % 可视化分段结果
    figure();
    
    % 颜色映射
    colors = lines(length(segments));
    
    % 绘制分段结果
    hold on;
    axis equal;
    for i = 1:length(segments)
        seg = segments(i);
        idx_range = seg.start:seg.end;
        
        if seg.type == 1
            % 长线段用黑色线绘制
            plot(x(idx_range), y(idx_range), 'LineWidth', 2, ...
                 'Color', 'black', 'DisplayName', '长直线段落');
            % 绘制拟合直线
            if isfield(seg, 'linear_params')
                x_fit = [x(seg.start), x(seg.end)];
                y_fit = seg.linear_params(1) * x_fit + seg.linear_params(2);
                plot(x_fit, y_fit, '--', 'Color', 'black', 'LineWidth', 1);
            end
            % 标记分段点
            plot(x(seg.start), y(seg.start), 'o', 'MarkerSize', 4, ...
                 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
        elseif seg.type == 2
            % 曲线段用普通线绘制
            plot(x(idx_range), y(idx_range), 'LineWidth', 2, ...
                 'Color', colors(i,:));
            % 标记分段点
            plot(x(seg.start), y(seg.start), 'o', 'MarkerSize', 4, ...
                 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
        elseif seg.type == 3
            % 转角用红色粗线绘制
            plot(x(idx_range), y(idx_range), 'LineWidth', 4, ...
                 'Color', 'r');
            % 标记分段点
            plot(x(seg.start), y(seg.start), 'o', 'MarkerSize', 4, ...
                 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
        end
    end
    title('分段结果');
    axis equal;
    grid on;
    hold off;

    % 输出分段信息
    fprintf('分段结果统计:\n');
    fprintf('总段数: %d\n', length(segments));
    for i = 1:length(segments)
        seg = segments(i);
        length_seg = seg.end - seg.start + 1;
        fprintf('段 %d: 点数=%d, 类型=%d\n', ...
                i, length_seg, seg.type);
    end
end