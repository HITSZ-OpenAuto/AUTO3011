close all;clc;clear;

% 曲面1：单位球 S(u,v)；法向与位置同向（可直接提供 Nfun）
S_circle = @(u,v) [cos(u).*sin(v); sin(u).*sin(v); cos(v)];
% Nfun = @(u,v) Sfun(u,v) ./ max(1e-15, vecnorm(Sfun(u,v))); % 可选
% uv-闭合曲线：纬线 v = const
Cuv_circle = @(t) [2*pi*t; pi/3 + 0*t];

% 曲面2：马鞍面
a = 1.2; b = 0.8;
S_hyper = @(u,v) [u; v; (u.^2)/(a^2) - (v.^2)/(b^2)];
% N_hyper = @(u,v) [ -2*u/(a^2);  2*v/(b^2); 1 ] ...
%                  ./ sqrt( (2*u/(a^2)).^2 + (2*v/(b^2)).^2 + 1 );
% uv-闭合曲线：uv-平面圆
R = 1.2;
Cuv_hyper = @(t) [ R*cos(2*pi*t); R*sin(2*pi*t) ];

% 曲面3：enneper 曲面
S_enneper = @(u,v) [u.*(1 - u.^2/3 + v.^2);...
                  v.*(1 - v^2/3 + u^2);...
                  u.^2 - v.^2];
% uv-闭合曲线：uv-平面圆，u, v ∈[-1.1,1.1]
R = 1.1;
Cuv_enneper = @(t) [ R*cos(2*pi*t); R*sin(2*pi*t) ];

% 指标
opts.Emax = 1e-3;  % 弦差
opts.Amax_deg = 10;  % 相邻线段夹角（角度制）
opts.Lmax = 0.2;  % 长度
opts.Bmax_deg = 2;  % 法向量夹角（角度制）
% opts.Nfun = Nfun;   % 有就用，没有就删掉这一行，自动差分

% 计时, 采样
tic
% opts.Name = "circle_1";
% [T,UV,P,N] = sample_surface_curve_normals(S_circle, Cuv_circle, opts);
% opts.Name = "hyper_1";
% [T,UV,P,N] = sample_surface_curve_normals(S_hyper, Cuv_hyper, opts);
opts.Name = "enneper_3";
[T,UV,P,N] = sample_surface_curve_normals(S_enneper, Cuv_enneper, opts);
toc
figure; plot3(P(1,:),P(2,:),P(3,:),'o-'); axis equal; grid on;

ax = gca;

% 绘制单位球上曲线小线段生成结果
% title_sphere = sprintf('Sphere + latitude curve with normals, N=%d', size(P,2));
% plot_surface_curve(S_circle, Cuv_circle, struct( ...
%     'tSamples', 800, ...
%     'showSurface', true, 'uRange', [0, 2*pi], 'vRange', [0, pi], ...
%     'nU', 80, 'nV', 40, 'normalEvery', 50, ...
%     'title', title_sphere,...
%     'ax', ax ));

% 绘制马鞍面上曲线小线段生成结果
% title_hyper = sprintf('Hyperbolic paraboloid + curve with normals, N=%d', size(P,2));
% plot_surface_curve(S_hyper, Cuv_hyper, struct( ...
%     'tSamples', 1200, ...
%     'showSurface', true, 'uRange', [-1.2,1.2], 'vRange', [-1.2,1.2], ...
%     'nU', 70, 'nV', 70, ...
%     'normalEvery', 50, 'normalScale', 0.25, ...
%     'title', title_hyper,...
%     'ax', ax));

% 绘制Enneper曲面上曲线小线段生成结果
title_enneper = sprintf('Enneper + curve with normals, N=%d', size(P,2));
plot_surface_curve(S_enneper, Cuv_enneper, struct( ...
    'tSamples', 1200, ...
    'showSurface', true, 'uRange', [-1.1,1.1], 'vRange', [-1.1,1.1], ...
    'nU', 70, 'nV', 70, ...
    'normalEvery', 50, 'normalScale', 0.25, ...
    'title', title_enneper,...
    'ax', ax));

function [T, UV, P, N] = sample_surface_curve_normals(Sfun, Cuv_fun, opts)
% 自适应最少点采样：
% 在曲面 S(u,v) 上，对 uv-平面里的闭合曲线映射得到的 3D 闭合曲线进行采样。
%
% 输入:
%   Sfun    : @(u,v)-> 3x1  曲面参数化
%   Cuv_fun : @(t)->  2x1  uv闭合曲线, t∈[0,1], Cuv_fun(0)=Cuv_fun(1)
%   opts:
%     必填: Emax, Amax_deg, Lmax, Bmax_deg, Name
%     可选: Nfun(u,v)->3x1（若未给则内部差分近似法向）
%           init_N=16, 初始采样步数
%           max_chord_check=30, 弦差检查采样数 
%           t_eps=1e-6, 最小时间间隔
%           du_rel=1e-6, dv_rel=1e-6, du_abs=1e-8, dv_abs=1e-8, 差分步
%
% 输出:
%   T : 1xN  采样的 t
%   UV: 2xN  对应 [u;v]
%   P : 3xN  对应 3D 点
%   N : 3xN  对应单位法向

    % ---------- 校验 & 默认 ----------
    assert(all(isfield(opts, {'Emax','Amax_deg','Lmax','Bmax_deg','Name'})), ...
        '需要提供 Emax, Amax_deg, Lmax, Bmax_deg');
    Amax = deg2rad(opts.Amax_deg);
    Bmax = deg2rad(opts.Bmax_deg);

    if ~isfield(opts,'init_N'),     opts.init_N = 16; end
    if ~isfield(opts,'max_chord_check'),  opts.max_chord_check = 30; end
    if ~isfield(opts,'max_iteration'),  opts.max_iteration = 20; end
    if ~isfield(opts,'max_trial'),  opts.max_trial = 10; end
    if ~isfield(opts,'t_eps'),      opts.t_eps = 1e-6; end

    if ~isfield(opts,'du_rel'), opts.du_rel = 1e-6; end
    if ~isfield(opts,'dv_rel'), opts.dv_rel = 1e-6; end
    if ~isfield(opts,'du_abs'), opts.du_abs = 1e-8; end
    if ~isfield(opts,'dv_abs'), opts.dv_abs = 1e-8; end

    haveNfun = isfield(opts,'Nfun') && ~isempty(opts.Nfun);
    if haveNfun
        Nfun = opts.Nfun;
    end

    % ---------- 便捷函数 ----------
    % 采样
    function [p,uv] = Cxyz(t)
        uv = Cuv_fun(t);
        p  = Sfun(uv(1), uv(2));
    end
    
    % 计算法向量
    function n = normal_uv(u,v, n_ref)
        if haveNfun
            n_raw = Nfun(u,v);
        else
            % 数值差分估计 Su, Sv
            du = max(opts.du_abs, opts.du_rel*(1+abs(u)));
            dv = max(opts.dv_abs, opts.dv_rel*(1+abs(v)));
            Su = (Sfun(u+du, v) - Sfun(u-du, v)) / (2*du);
            Sv = (Sfun(u, v+dv) - Sfun(u, v-dv)) / (2*dv);
            n_raw = cross(Su, Sv);
            if norm(n_raw) < 1e-14
                % 退一步：放大步长再尝试
                Su = (Sfun(u+10*du, v) - Sfun(u-10*du, v)) / (20*du);
                Sv = (Sfun(u, v+10*dv) - Sfun(u, v-10*dv)) / (20*dv);
                n_raw = cross(Su, Sv);
            end
        end
        if norm(n_raw) < 1e-14
            % 退化情形：回退到参考法向或默认 Z 轴
            if nargin>=3 && ~isempty(n_ref), n = n_ref;
            else, n = [0;0;1];
            end
            return;
        end
        n = n_raw / norm(n_raw);
        % 与参考法向做方向连续化（防止 ±n 翻转）
        if nargin>=3 && ~isempty(n_ref)
            if dot(n, n_ref) < 0, n = -n; end
        end
    end

    function [p,uv,n] = sample_at_t(t, n_ref)
        [p,uv] = Cxyz(t);
        n = normal_uv(uv(1), uv(2), n_ref);
    end

    function ang = safe_angle(a,b)
        na = norm(a); nb = norm(b);
        if na==0 || nb==0, ang = 0; return; end
        c = dot(a,b)/(na*nb);
        c = max(-1,min(1,c));
        ang = acos(c);
    end

    function d = point_to_seg_dist(P, A, B)
        AB = B - A;  ab2 = dot(AB,AB);
        if ab2 < eps, d = norm(P-A); return; end
        location = dot(P - A, AB) / ab2;
        location = max(0,min(1,location));
        Q = A + location*AB; % 投影点（垂足）
        d = norm(P - Q);
    end

    function [bad, maxd] = chord_exceed(ta, tb, pa, pb)
        tm = linspace(ta, tb, opts.max_chord_check);
        maxd = 0;
        % 采样计算最大偏差距离
        for i=2:length(tm)-1
            t_sample=tm(i);
            pm = Cxyz(t_sample); pm = pm(:,1); % 只取点
            d = point_to_seg_dist(pm, pa, pb);
            if maxd < d, maxd = d; end
            if d > opts.Emax, bad = true; return; end
        end
        bad = false;
    end

    % 重点: 检查每一段是否符合要求
    function [ok, metric] = seg_ok(ta, tb, prev_vec, prev_n, print_log, metric)
        [pa,~,na] = sample_at_t(ta, prev_n);    % 端点A，法向与 prev_n 对齐
        [pb,~,nb] = sample_at_t(tb, na);        % 端点B，法向与 na 连续

        ok = true;
        % 长度约束
        seg_len = norm(pb - pa);
        if seg_len > opts.Lmax + 1e-14
            ok = false; 
            if nargin >=5 && print_log
                fprintf("seg_len = %f exceed limit, " + ...
                    "ta = %f, tb = %f\n", seg_len, ta, tb);
            end
        end

        % 相邻线段夹角约束
        if ~isempty(prev_vec)
            ang = safe_angle(prev_vec, pb - pa);
            if ang > Amax + 1e-14
                ok = false; 
                if nargin >=5 && print_log
                    fprintf("ang = %f exceed limit, " + ...
                        "ta = %f, tb = %f\n", ang, ta, tb);
                end
            end
        end

        % 法向夹角约束
        bang = safe_angle(na, nb);
        if bang > Bmax + 1e-14
            ok = false; 
            if nargin >=5 && print_log
                fprintf("bang %f exceed limit, " + ...
                        "ta = %f, tb = %f\n", bang, ta, tb);
            end
        end

        % 弦差（曲线 vs 直线段）
        [bad, chord]=chord_exceed(ta, tb, pa, pb);
        if bad
            ok = false; 
            if nargin >=5 && print_log
                fprintf("chord exceed limit, " + ...
                        "ta = %f, tb = %f\n", ta, tb);
            end
        end

        % 保存指标
        if nargin == 6
            new_metric = struct('chord', chord, 'angle', rad2deg(ang), ...
                'length', seg_len, 'angle_of_norm', rad2deg(bang));
            metric = [metric, new_metric];
        end
    end

    function t1 = binary_max_feasible(t0_, low_, high_, prev_vec_, prev_n_)
        L = low_; H = high_;
        for k = 1:opts.max_iteration
            if (H - L) < opts.t_eps, break; end
            mid = 0.5*(L + H);
            if seg_ok(t0_, mid, prev_vec_, prev_n_)
                L = mid; 
            else
                H = mid; 
            end
        end
        t1 = L;
        if ~(seg_ok(t0_, L, prev_vec_, prev_n_))
            error("[binary_max_feasible] 选取的时间不合理！");
        end
    end

    % 在首点处保证：段角 ≤ Amax、法向角 ≤ Bmax
    function [T,UV,P,N] = fix_closure(T,UV,P,N)
        max_iter = 60; it = 0;
        while it < max_iter
            it = it + 1;
            if size(P,2) < 3, break; end
            v_last  = P(:,1) - P(:,end);  % 第一个点及其之前
            v_first = P(:,2) - P(:,1);    % 第一个点及其之后
            ang_v = safe_angle(v_last, v_first);
            ang_n = safe_angle(N(:,end), N(:,1));

            if ang_v <= Amax + 1e-14 && ang_n <= Bmax + 1e-14
                break; % 满足闭合约束
            end

            % 选择需要细分的一侧（优先改进超标更多的一项）
            if ang_n - Bmax > ang_v - Amax
                % 优先修法向夹角：在法向变化更大的那一段细分
                % 比较 N(end-1)->N(end) 与 N(1)->N(2) 的差异
                ang_lastN  = safe_angle(N(:,end-1), N(:,end));
                ang_firstN = safe_angle(N(:,1),     N(:,2));
                split_last = (ang_lastN >= ang_firstN);
            else
                % 优先修折线夹角：细分更长的那段
                split_last = (norm(v_last) >= norm(v_first));
            end

            if split_last
                t_mid = 0.5*(T(end-1) + T(end));
                [p_mid, uv_mid, n_mid] = sample_at_t(t_mid, N(:,end-1));
                T  = [T(1:end-1), t_mid, T(end)];
                UV = [UV(:,1:end-1), uv_mid, UV(:,end)];
                P  = [P(:,1:end-1),  p_mid,  P(:,end)];
                N  = [N(:,1:end-1),  n_mid,  N(:,end)];
            else
                t_mid = 0.5*(T(1) + T(2));
                [p_mid, uv_mid, n_mid] = sample_at_t(t_mid, N(:,1));
                T  = [T(1), t_mid, T(2:end)];
                UV = [UV(:,1), uv_mid, UV(:,2:end)];
                P  = [P(:,1),  p_mid,  P(:,2:end)];
                N  = [N(:,1),  n_mid,  N(:,2:end)];
            end
        end
    end

    % 进行总体检验与数据保存
    function ok = post_checkout(T,N,P)
        check_prev_n = N(:,end); check_prev_vec = P(:,1) - P(:,end);
        metric = [];
        for i = 1:length(T) - 1
            [ok, metric] = seg_ok(T(i), T(i+1), check_prev_vec, check_prev_n, true, metric);
            [pb,~,check_prev_n] = sample_at_t(T(i+1), check_prev_n);
            [pa,~,~] = sample_at_t(T(i), check_prev_n);
            check_prev_vec = pb - pa;
        end
        metric_table = struct2table(metric);
        format long
        max(metric_table{:,1:4})
        % 保存数据供数据分析
        writetable(metric_table, "metric_" + opts.Name + ".csv");
    end

    % ---------- 初始化 ----------
    t0 = 0.0;
    [p0,uv0,n0] = sample_at_t(t0, []);
    T  = t0;   P  = p0;   UV = uv0;   N = n0;

    prev_vec    = [];     % 上一段方向
    prev_n      = n0;     % 上一点法向
    base_dt     = 1/opts.init_N;
    t           = t0;
    trial       = 0;

    % ---------- 主循环 ----------
    while t < 1 - opts.t_eps
        dt_try = min(base_dt, 1 - t);
        ok = seg_ok(t, t+dt_try, prev_vec, prev_n);
        if ~ok
            shrink = 0;
            while ~ok && dt_try > opts.t_eps && shrink < opts.max_iteration
                dt_try = dt_try * 0.5;
                ok = seg_ok(t, t+dt_try, prev_vec, prev_n);
                shrink = shrink + 1;
            end
            if ~ok
                fprintf('无法从 t=%.6f 前进。尝试调整……\n', t);
                trial = trial + 1;
                if trial == 1
                    gap = t - T(end-1);
                    % 清除原来的末点
                    T(:,end) = []; P(:,end) = []; 
                    UV(:,end) = []; N(:,end) = [];
                end
                if trial >= opts.max_trial
                    if isempty(T)
                        error('约束过严，无法前进。');
                    else
                        fprintf('继续调整上一个点！'); % 尚未经过充分测试
                        gap = t - T(end-1);
                        T(:,end) = []; P(:,end) = []; 
                        UV(:,end) = []; N(:,end) = [];
                        trial = 0;
                    end
                end
                t = t - gap/opts.max_trial;
                [p2,~,n2] = sample_at_t(t, N(:,end));
                prev_vec = p2 - P(:,end);
                prev_n = n2;
                continue
            else
                if trial > 0
                    fprintf('调整成功！现在t=%.6f \n', t);
                    trial = 0;
                    [p2,uv2,n2] = sample_at_t(t, N(:,end));
                    % 将调整点加入
                    T  = [T,  t];
                    P  = [P,  p2];
                    UV = [UV, uv2];
                    N  = [N,  n2];
                    continue
                end
            end
        end

        last_good = t + dt_try;
        dt_exp    = 0.5 * dt_try;
        first     = true;
        % 指数扩张，直到首次不可行
        while true
            dt_exp = min(2*dt_exp, 1 - t);
            if dt_exp <= 0, break; end
            cand = t + dt_exp;
            ok = seg_ok(t, cand, prev_vec, prev_n);
            if ok
                last_good = cand;
                first = false;
                if abs(cand - 1) < opts.t_eps, break; end
            else
                if first
                    error("[exp_expand] last_good 不合理！");
                end
                break;
            end
            if abs(dt_exp - (1 - t)) < opts.t_eps, break; end
        end

        % 用二分逼近最大可行 t1
        high = min(1, t + dt_exp);
        t1 = binary_max_feasible(t, last_good, high, prev_vec, prev_n);
%         fprintf("t = %f, last_good = %f, high = %f, t1 = %f\n", t, last_good, high, t1);

        % 接受 t1
        [p1,uv1,n1] = sample_at_t(t1, prev_n);
        T  = [T,  t1];
        P  = [P,  p1];
        UV = [UV, uv1];
        N  = [N,  n1];

        prev_vec = p1 - P(:,end-1);
        prev_n   = n1;
        t = t1;
    end
    
    % ---------- 闭合修正 ----------
    [T,UV,P,N] = fix_closure(T,UV,P,N);

    % 轨迹检查及数据分析
    if post_checkout(T, N, P)
        disp("检查通过！");
        % 保存结果
        results = table(P(1,:)', P(2,:)', P(3,:)', N(1,:)', N(2,:)', N(3,:)');
        writetable(results, "result_" + opts.Name + ".txt", ...
            'WriteVariableNames',false,'Delimiter',' ');
    end
end

function [ax, out] = plot_surface_curve(Sfun, Cuv_fun, opts)
% 根据 Sfun(u,v) 与 Cuv_fun(t) 绘制曲面上的闭合/开曲线，并可选绘制曲面与法向
%
% 输入
%   Sfun    : @(u,v)->3x1  曲面参数化（返回 [x;y;z]）
%   Cuv_fun : @(t)->2x1    uv 参数平面中的曲线，t ∈ [0,1]
%   opts    : 结构体，可选字段（都有默认值）：
%       % 曲线采样
%       .tSamples = 1000;           % 曲线均匀采样点数
%       .close    = true;           % 若首尾不重合，是否重复首点闭合
%       .curveLineWidth = 2;
%       .curveColor     = [0.85 0.2 0.1];
%
%       % 曲面网格（可选）
%       .showSurface = false;       % 是否绘制曲面
%       .uRange = [-1, 1];          % u 范围
%       .vRange = [-1, 1];          % v 范围
%       .nU = 60; .nV = 60;         % 曲面网格密度
%       .surfFaceAlpha = 0.15;
%       .surfEdgeAlpha = 0.2;
%       .surfEdgeColor = [0.5 0.5 0.5];
%
%       % 法向（可选）
%       .Nfun = [];                 % @(u,v)->3x1 单位法向；不给则不画法向
%       .normalEvery  = -1;         % 每隔多少个点画一次法向, -1 表示不启用
%       .normalScale  = 0.2;        % 法向箭头长度缩放
%       .normalColor  = [0.1 0.3 0.9];
%
%       % 轴/窗口
%       .ax = [];                   % 传入已有坐标轴；为空则新建
%       .title = 'Surface curve';
%
% 输出
%   ax  : 坐标轴句柄
%   out : 结构体，便于后续使用
%       .t  (1xN)     均匀采样的参数
%       .UV (2xN)     曲线对应的 (u,v)
%       .P  (3xN)     空间曲线点
%       .h  (struct)  图元句柄（curve/surface/normals）

    % ---------- 默认参数 ----------
    if nargin < 3, opts = struct(); end
    setdef = @(s, f, v) ( isfield(s,f) && ~isempty(s.(f)) ) || ~isempty(setfield(s,f,[]));
    if ~isfield(opts,'tSamples'),      opts.tSamples = 1000; end
    if ~isfield(opts,'close'),         opts.close    = true; end
    if ~isfield(opts,'curveLineWidth'),opts.curveLineWidth = 2; end
    if ~isfield(opts,'curveColor'),    opts.curveColor     = [0.85 0.2 0.1]; end

    if ~isfield(opts,'showSurface'),   opts.showSurface = false; end
    if ~isfield(opts,'uRange'),        opts.uRange = [-1, 1]; end
    if ~isfield(opts,'vRange'),        opts.vRange = [-1, 1]; end
    if ~isfield(opts,'nU'),            opts.nU = 60; end
    if ~isfield(opts,'nV'),            opts.nV = 60; end
    if ~isfield(opts,'surfFaceAlpha'), opts.surfFaceAlpha = 0.15; end
    if ~isfield(opts,'surfEdgeAlpha'), opts.surfEdgeAlpha = 0.2; end
    if ~isfield(opts,'surfEdgeColor'), opts.surfEdgeColor = [0.5 0.5 0.5]; end

    if ~isfield(opts,'Nfun'),          opts.Nfun = []; end
    if ~isfield(opts,'normalEvery'),   opts.normalEvery = -1; end
    if ~isfield(opts,'normalScale'),   opts.normalScale = 0.2; end
    if ~isfield(opts,'normalColor'),   opts.normalColor = [0.1 0.3 0.9]; end

    if ~isfield(opts,'ax') || isempty(opts.ax)
        fig = figure('Color','w'); ax = axes('Parent',fig); 
    else
        ax = opts.ax;
    end
    if ~isfield(opts,'title'), opts.title = 'Surface curve'; end

    hold(ax,'on'); grid(ax,'on');

    % ---------- 曲线均匀采样 ----------
    t  = linspace(0,1,opts.tSamples);
    UV = zeros(2, numel(t));
    P  = zeros(3, numel(t));
    for i = 1:numel(t)
        uv = Cuv_fun(t(i));
        UV(:,i) = uv(:);
        p = Sfun(uv(1), uv(2));
        P(:,i) = p(:);
    end

    % 闭合处理（可选）：末点重复首点，便于显示闭合折线
    if opts.close
        if norm(P(:,end) - P(:,1)) > 1e-12
            t  = [t, 1];
            UV = [UV, UV(:,1)];
            P  = [P,  P(:,1)];
        else
            t(end) = 1; % 保证参数末端对齐
        end
    end

    % ---------- 曲面（可选） ----------
    hSurf = [];
    if opts.showSurface
        u_lin = linspace(opts.uRange(1), opts.uRange(2), opts.nU);
        v_lin = linspace(opts.vRange(1), opts.vRange(2), opts.nV);
        [U,V] = meshgrid(u_lin, v_lin);
        X = zeros(size(U)); Y = X; Z = X;
        for i = 1:numel(U)
            p = Sfun(U(i), V(i));
            X(i) = p(1); Y(i) = p(2); Z(i) = p(3);
        end
        hSurf = surf(ax, X, Y, Z, 'FaceAlpha', opts.surfFaceAlpha, ...
            'EdgeAlpha', opts.surfEdgeAlpha, 'EdgeColor', opts.surfEdgeColor, ...
            'FaceColor',[0.8 0.8 0.85]);
    end

    % ---------- 曲线 ----------
    hCurve = plot3(ax, P(1,:), P(2,:), P(3,:), '-', ...
        'LineWidth', opts.curveLineWidth, 'Color', opts.curveColor);

    % ---------- 法向（可选） ----------
    hQ = [];
    if opts.normalEvery > 0
        idx = 1:opts.normalEvery:numel(t);
        Q  = zeros(3, numel(idx));
        for k = 1:numel(idx)
            i = idx(k);
            if ~isempty(opts.Nfun)
                n = opts.Nfun(UV(1,i), UV(2,i));
            else
                du_rel = 1e-6; dv_rel = 1e-6; du_abs = 1e-8; dv_abs = 1e-8; 
                u = UV(1,i);
                v = UV(2,i);
                du = max(du_abs, du_rel*(1+abs(u)));
                dv = max(dv_abs, dv_rel*(1+abs(v)));
                Su = (Sfun(u+du, v) - Sfun(u-du, v)) / (2*du);
                Sv = (Sfun(u, v+dv) - Sfun(u, v-dv)) / (2*dv);
                n = cross(Su, Sv);
            end
            if ~isempty(n)
                n = n(:) / max(1e-15, norm(n));
            else
                n = [0;0;1];
            end
            Q(:,k) = n;
        end
        scale = opts.normalScale;
        hQ = quiver3(ax, P(1,idx), P(2,idx), P(3,idx), ...
                          Q(1,:),   Q(2,:),   Q(3,:), scale, ...
                          'Color', opts.normalColor, 'LineWidth', 1.2);        
    end

    % ---------- 轴样式 ----------
    axis(ax,'equal'); view(ax, 45, 25);
    xlabel(ax,'X'); ylabel(ax,'Y'); zlabel(ax,'Z');
    title(ax, opts.title);

    % ---------- 输出 ----------
    out.t  = t;
    out.UV = UV;
    out.P  = P;
    out.h.surface = hSurf;
    out.h.curve   = hCurve;
    out.h.normals = hQ;
end
