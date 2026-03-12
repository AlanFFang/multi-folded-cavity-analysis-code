%% 三折谐振腔优化程序
clear; clc;

% --- 1. 基础物理参数 ---
r1 = 30e-3; r4 = 300e-3; d = 10e-3; 
f_target = 20e6; omega0 = 2*pi*f_target;
sigma = 5.8e7; mu0 = 4*pi*1e-7; epi = 8.854e-12; 
eta = sqrt(mu0/epi); c = 3e8; l = c/f_target/12;
Rs_surf_ref = sqrt(omega0*mu0/(2*sigma));

% --- 2. 扫描网格 ---
n_grid = 100; 
r2_list = linspace(r1+d+0.002, r4-2*d-0.01, n_grid);
r3_list = linspace(r1+2*d+0.01, r4-d-0.002, n_grid);
[R2, R3] = meshgrid(r2_list, r3_list);

Rs_matrix = nan(size(R2));
Freq_matrix = nan(size(R2));

% --- 3. 核心计算循环 ---
for i = 1:size(R2,1)
    for j = 1:size(R2,2)
        curr_r2 = R2(i,j); curr_r3 = R3(i,j);
        if curr_r3 > (curr_r2 + d + 0.001)
            % 定义寻根目标：Im(1/Zin) = 0
            target_f = @(w) get_ImY(w, curr_r2, curr_r3, r1, r4, d, sigma, mu0, eta, c, l);
            try
                w_res = fzero(target_f, 2*pi*f_target, optimset('Display','off'));
                [~, Zin] = get_ImY(w_res, curr_r2, curr_r3, r1, r4, d, sigma, mu0, eta, c, l);
                Rs_matrix(i,j) = real(Zin);
                Freq_matrix(i,j) = w_res / (2*pi);
            catch
                continue;
            end
        end
    end
end

% --- 4. 寻找 20MHz 附近的最高阻抗解 ---
f_tol = 0.05e6; 
valid_idx = find(abs(Freq_matrix - f_target) < f_tol); % 定义 valid_idx

% --- 5. 绘图与平滑处理 ---

%% ================== 1. 核心计算：寻找最佳工程点 ==================
% 提取 20MHz 轨迹并执行 [r2:30-100, r3:130-290] 区间筛选
C = contourc(r2_list, r3_list, Freq_matrix/1e6, [20 20]);
raw_x = C(1, 2:end); raw_y = C(2, 2:end);
mask = (raw_x >= 30e-3) & (raw_x <= 100e-3) & (raw_y >= 130e-3) & (raw_y <= 290e-3);
line_x = raw_x(mask); line_y = raw_y(mask);

% 计算 Rs 和 分量敏感度 (导纳斜率法)
Rs_line = interp2(R2, R3, Rs_matrix/1000, line_x, line_y, 'linear');
df_dr2 = zeros(size(line_x)); df_dr3 = zeros(size(line_x));
dr = 0.01e-3; df_step = 100;

for k = 1:length(line_x)
    cx = line_x(k); cy = line_y(k);
    dB_df = (get_ImY(2*pi*(20e6+df_step), cx, cy, r1, r4, d, sigma, mu0, eta, c, l) - ...
             get_ImY(2*pi*(20e6-df_step), cx, cy, r1, r4, d, sigma, mu0, eta, c, l)) / (2*df_step);
    dB_dr2 = (get_ImY(2*pi*20e6, cx+dr, cy, r1, r4, d, sigma, mu0, eta, c, l) - ...
              get_ImY(2*pi*20e6, cx-dr, cy, r1, r4, d, sigma, mu0, eta, c, l)) / (2*dr*1000);
    dB_dr3 = (get_ImY(2*pi*20e6, cx, cy+dr, r1, r4, d, sigma, mu0, eta, c, l) - ...
              get_ImY(2*pi*20e6, cx, cy-dr, r1, r4, d, sigma, mu0, eta, c, l)) / (2*dr*1000);
    df_dr2(k) = -(dB_dr2 / dB_df) / 1000; 
    df_dr3(k) = -(dB_dr3 / dB_df) / 1000;
end

% 综合评价函数选取最佳点 (Weight: Rs 0.6, Sensitivity 0.4)
total_S = sqrt(df_dr2.^2 + df_dr3.^2);
norm_Rs = (Rs_line - min(Rs_line)) / (max(Rs_line) - min(Rs_line) + 1e-6);
norm_S = (total_S - min(total_S)) / (max(total_S) - min(total_S) + 1e-6);
score = 0.6 * norm_Rs - 0.4 * norm_S;
[~, best_idx] = max(score);

%% ================== 2. 绘图 1：全局细节云图 ==================
[R2_f, R3_f] = meshgrid(linspace(min(r2_list),max(r2_list),400), linspace(min(r3_list),max(r3_list),400));
Rs_smooth = interp2(R2, R3, Rs_matrix, R2_f, R3_f, 'linear');
Freq_smooth = interp2(R2, R3, Freq_matrix, R2_f, R3_f, 'linear');

figure('Color','w','Position',[100 100 900 700]);
hold on;
set(gca, 'FontSize', 14, 'LineWidth', 1.2); 
% A. 背景 Rs 云图
[~, h_bg] = contourf(R2_f*1000, R3_f*1000, 2*Rs_smooth/1000, 50, 'LineStyle', 'none');
cb = colorbar; cb.Label.String = 'Shunt Impedance R_s (k\Omega)';
cb.Label.FontSize = 15; % 加大 Colorbar 标签字体

% B. 辅助频率线 (白色虚线)
f_lines = [15, 17.5, 22.5, 25];
[C_ref, h_ref] = contour(R2_f*1000, R3_f*1000, Freq_smooth/1e6, f_lines, 'w--', 'LineWidth', 1.0);
ht1 = clabel(C_ref, h_ref, 'FontSize', 12, 'Color', 'w', 'FontWeight', 'bold'); % 加大等高线数值
for k=1:length(ht1), set(ht1(k), 'String', [get(ht1(k),'String'), ' MHz']); end

% C. 目标 20MHz 线 (加厚黑实线 + 红字白底标签)
[C_t, h_t] = contour(R2_f*1000, R3_f*1000, Freq_smooth/1e6, [20 20], 'k-', 'LineWidth', 3);
ht2 = clabel(C_t, h_t, 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold'); % 重点突出 20MHz 标签
for k=1:length(ht2)
    set(ht2(k), 'String', '20 MHz', 'BackgroundColor', 'w', 'Margin', 2, 'EdgeColor', 'r');
end

% D. 仅标注最佳工程点 (青色方块，带外框)
plot(line_x(best_idx)*1000, line_y(best_idx)*1000, 'rp', ...
     'MarkerSize', 20, 'LineWidth', 1.5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y');

% E. 装饰
xlabel('r_2 (mm)', 'FontSize', 16, 'FontWeight', 'bold'); 
ylabel('r_3 (mm)', 'FontSize', 16, 'FontWeight', 'bold');
title('Triple-folded Cavity Parameter Scaning', 'FontSize', 18);
grid on; colormap parula;
hold off;

%% ================== 3. 绘图 2：局部平衡切片图 ==================
figure('Color','w','Position',[100 100 900 700]);
ax = gca; hold on;
set(gca, 'FontSize', 14, 'LineWidth', 1.2); 

yyaxis left
h1 = plot(line_x*1000, 2*Rs_line, 'k-', 'LineWidth', 3, 'DisplayName', 'R_s');
% 标注最佳点在 Rs 曲线上的位置
plot(line_x(best_idx)*1000, 2*Rs_line(best_idx), 'kp', 'MarkerSize', 16, 'MarkerFaceColor', 'y', 'HandleVisibility', 'off');
ylabel('Shunt Impedance R_s (k\Omega)', 'FontSize', 16, 'FontWeight', 'bold');
ax.YColor = [0 0 0];

yyaxis right
h2 = plot(line_x*1000, df_dr2, 'r--', 'LineWidth', 2, 'DisplayName', '\partialf / \partialr_2');
h3 = plot(line_x*1000, df_dr3, 'b--', 'LineWidth', 2, 'DisplayName', '\partialf / \partialr_3');
% 标注分量敏感度
plot(line_x(best_idx)*1000, df_dr2(best_idx), 'ro', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
plot(line_x(best_idx)*1000, df_dr3(best_idx), 'bo', 'MarkerFaceColor', 'b', 'HandleVisibility', 'off');

ylabel('Frequency Sensitivity (kHz/mm)', 'FontSize', 16, 'FontWeight', 'bold');
yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');

grid on;
xlabel('Radius r_2 (mm) along 20 MHz trajectory', 'FontSize', 16);
title('Shunt Impedance vs. Frequency Stability in 20MHz contour line', 'FontSize', 18);
legend([h1, h2, h3], 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 16);

%% ================== 4. 最佳工程点全参数核算 ==================

%% ================== 4. 最佳工程点全参数核算 (数据对齐版) ==================

% 1. 从插值数据中提取 Rs (确保和图二标注一致)
Rs_from_plot = Rs_line(best_idx); % 单位: kOhm

% 2. 物理校准计算 (为了得到 Q 和 R/Q)
r2_best = line_x(best_idx);
r3_best = line_y(best_idx);

% 寻找该尺寸下的精确谐振频率 (确保虚部归零)
f_res = fzero(@(f) get_ImY(2*pi*f, r2_best, r3_best, r1, r4, d, sigma, mu0, eta, c, l), 20e6);

% 在谐振频率下计算 Q 值
dw = 10;
[~, Zin_res] = get_ImY(2*pi*f_res, r2_best, r3_best, r1, r4, d, sigma, mu0, eta, c, l);
ImY_up = get_ImY(2*pi*f_res + dw, r2_best, r3_best, r1, r4, d, sigma, mu0, eta, c, l);
ImY_dn = get_ImY(2*pi*f_res - dw, r2_best, r3_best, r1, r4, d, sigma, mu0, eta, c, l);
dB_dw = (ImY_up - ImY_dn) / (2 * dw);
G_res = real(1/Zin_res);
Q_val = (2*pi*f_res / 2) * dB_dw / G_res;

% 计算理论 Rs (kOhm) 用于对比
Rs_theoretical = real(Zin_res) / 1000;

% 计算几何阻抗 R/Q (Ohm)
RoQ_val = (Rs_from_plot * 1000) / Q_val; % 这里用图上的 Rs 来算 R/Q 保证一致性

% 3. 打印对齐后的结果
fprintf('\n================== 最佳工程点参数 (对齐版) ==================\n');
fprintf('几何尺寸:\n');
fprintf('  r2 = %.2f mm, r3 = %.2f mm\n', r2_best*1000, r3_best*1000);
fprintf('\n电磁参数 (用于海报展示):\n');
fprintf('  谐振频率 f_res:    %.4f MHz\n', f_res/1e6);
fprintf('  分流阻抗 Rs:       %.2f kOhm \n', Rs_from_plot);
fprintf('  品质因数 Q:        %.1f\n', Q_val);
fprintf('  几何阻抗 R/Q:      %.2f Ohm\n', RoQ_val);
fprintf('\n数据一致性检查:\n');
fprintf('  图表 Rs: %.2f kOhm vs. 理论 Rs: %.2f kOhm\n', Rs_from_plot, Rs_theoretical);
fprintf('============================================================\n');
% [R2_f, R3_f] = meshgrid(linspace(min(r2_list),max(r2_list),400), linspace(min(r3_list),max(r3_list),400));
% Rs_smooth = interp2(R2, R3, Rs_matrix, R2_f, R3_f, 'linear');
% Freq_smooth = interp2(R2, R3, Freq_matrix, R2_f, R3_f, 'linear');
% 
% figure('Color','w','Position',[100 100 900 700]);
% hold on;
% [~, h_bg] = contourf(R2_f*1000, R3_f*1000, Rs_smooth/1000, 50, 'LineStyle', 'none');
% 
% % 辅助频率线标注
% f_lines = [15, 17.5, 22.5, 25];
% [C_ref, h_ref] = contour(R2_f*1000, R3_f*1000, Freq_smooth/1e6, f_lines, 'w--', 'LineWidth', 1.0);
% ht1 = clabel(C_ref, h_ref, 'FontSize', 9, 'Color', 'w');
% for k=1:length(ht1), set(ht1(k), 'String', [get(ht1(k),'String'), ' MHz']); end
% 
% % 目标 20MHz 线标注 (红字白底)
% [C_t, h_t] = contour(R2_f*1000, R3_f*1000, Freq_smooth/1e6, [20 20], 'k-', 'LineWidth', 3);
% ht2 = clabel(C_t, h_t, 'FontSize', 12, 'Color', 'r', 'FontWeight', 'bold');
% for k=1:length(ht2)
%     set(ht2(k), 'String', '20 MHz', 'BackgroundColor', 'w', 'Margin', 2, 'EdgeColor', 'r');
% end
% 
% % 绘制最优解点
% if ~isempty(valid_idx)
%     [max_R, local_idx] = max(Rs_matrix(valid_idx));
%     best_ptr = valid_idx(local_idx);
%     plot(R2(best_ptr)*1000, R3(best_ptr)*1000, 'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'y');
% end
% 
% xlabel('r_2 (mm)'); ylabel('r_3 (mm)'); title('Resonant Frequency and Shunt Impedance Optimization');
% colorbar; colormap parula; grid on; hold off;
% % --- 6. 打印最优参数结果 ---
% if ~isempty(valid_idx)
%     % 已经在之前的逻辑中找到了 best_ptr
%     best_r2_mm = R2(best_ptr) * 1000;
%     best_r3_mm = R3(best_ptr) * 1000;
%     best_Rs_kOhm = Rs_matrix(best_ptr) / 1000;
%     best_Freq_MHz = Freq_matrix(best_ptr) / 1e6;
% 
%     fprintf('\n================ 最优设计参数 (20MHz 约束下) ================ \n');
%     fprintf('最大分流阻抗 (Rs):      %.2f kOhm\n', best_Rs_kOhm);
%     fprintf('最优中间半径 (r2):      %.2f mm\n',   best_r2_mm);
%     fprintf('最优中间半径 (r3):      %.2f mm\n',   best_r3_mm);
%     fprintf('精确谐振频率 (f_res):   %.4f MHz\n',  best_Freq_MHz);
%     fprintf('============================================================ \n');
% else
%     fprintf('\n[警告] 在设定的频率容差 (f_tol) 范围内未找到谐振点，请扩大扫描范围或减小 f_tol。\n');
% end
% 
% 
% %% --- 9. 性能与稳定性全参数对比 (Rs, df/dr2, df/dr3) ---
% 
% % 1. 提取 20MHz 轨迹点
% C = contourc(r2_list, r3_list, Freq_matrix/1e6, [20 20]);
% raw_x = C(1, 2:end); raw_y = C(2, 2:end);
% mask = (raw_x >= 30e-3) & (raw_x <= 100e-3) & (raw_y >= 130e-3) & (raw_y <= 290e-3);
% line_x = raw_x(mask); line_y = raw_y(mask);
% 
% % 2. 预分配并计算
% Rs_line = interp2(R2, R3, Rs_matrix/1000, line_x, line_y, 'linear');
% df_dr2 = zeros(size(line_x));
% df_dr3 = zeros(size(line_x));
% dr = 0.01e-3; df = 100;
% 
% for k = 1:length(line_x)
%     cx = line_x(k); cy = line_y(k);
%     % 频率斜率 (用于归一化)
%     ImY_p = get_ImY(2*pi*(20e6+df), cx, cy, r1, r4, d, sigma, mu0, eta, c, l);
%     ImY_m = get_ImY(2*pi*(20e6-df), cx, cy, r1, r4, d, sigma, mu0, eta, c, l);
%     dB_df = (ImY_p - ImY_m) / (2*df);
% 
%     % r2, r3 扰动
%     dB_dr2 = (get_ImY(2*pi*20e6, cx+dr, cy, r1, r4, d, sigma, mu0, eta, c, l) - ...
%               get_ImY(2*pi*20e6, cx-dr, cy, r1, r4, d, sigma, mu0, eta, c, l)) / (2*dr*1000);
%     dB_dr3 = (get_ImY(2*pi*20e6, cx, cy+dr, r1, r4, d, sigma, mu0, eta, c, l) - ...
%               get_ImY(2*pi*20e6, cx, cy-dr, r1, r4, d, sigma, mu0, eta, c, l)) / (2*dr*1000);
% 
%     df_dr2(k) = -(dB_dr2 / dB_df) / 1000; % kHz/mm
%     df_dr3(k) = -(dB_dr3 / dB_df) / 1000;
% end
% 
% % --- 3. 绘图：双 Y 轴对比 ---
% figure('Color','w','Position',[100 100 900 600]);
% ax = gca; hold on;
% 
% % 左轴：分流阻抗 (Rs)
% yyaxis left
% h1 = plot(line_x*1000, Rs_line, 'k-', 'LineWidth', 3, 'DisplayName', 'Shunt Impedance R_s');
% ylabel('Shunt Impedance R_s (k\Omega)', 'FontSize', 14, 'FontWeight', 'bold');
% ax.YColor = [0 0 0]; % 黑色
% 
% % 右轴：敏感度分量 (df/dr)
% yyaxis right
% h2 = plot(line_x*1000, df_dr2, 'r--', 'LineWidth', 2, 'DisplayName', '\partialf / \partialr_2');
% h3 = plot(line_x*1000, df_dr3, 'b--', 'LineWidth', 2, 'DisplayName', '\partialf / \partialr_3');
% ylabel('Sensitivity (kHz/mm)', 'FontSize', 14, 'FontWeight', 'bold');
% yline(0, 'k:', 'HandleVisibility', 'off'); % 零线
% 
% % 格式美化
% grid on;
% xlabel('Radius r_2 (mm) along 20 MHz ridge', 'FontSize', 14);
% title('Stability-Efficiency Trade-off Analysis', 'FontSize', 16);
% legend([h1, h2, h3], 'Location', 'southoutside', 'Orientation', 'horizontal');
% 
% % --- 4. 寻找"工程黄金点" ---
% % 在 Rs 高于均值的点中，寻找总敏感度模值最小的点
% total_S = sqrt(df_dr2.^2 + df_dr3.^2);
% [min_S, best_idx] = min(total_S);
% 
% plot(line_x(best_idx)*1000, df_dr2(best_idx), 'ro', 'MarkerFaceColor', 'r');
% plot(line_x(best_idx)*1000, df_dr3(best_idx), 'bo', 'MarkerFaceColor', 'b');
% fprintf('\n--- 建议工程设计点 ---\n');
% fprintf('r2 = %.2f mm, r3 = %.2f mm\n', line_x(best_idx)*1000, line_y(best_idx)*1000);
% fprintf('Rs = %.2f kOhm\n', Rs_line(best_idx));
% fprintf('Total Sensitivity = %.2f kHz/mm\n', min_S);


% --- 辅助函数 ---
function [ImY, Zin] = get_ImY(w, r2, r3, r1, r4, d, sigma, mu0, eta, c, l)
    Rs = sqrt(w * mu0 / (2 * sigma));
    % 特征阻抗及损耗
    ln3 = log((r2-d/2)/(r1+d/2)); ln2 = log((r3-d/2)/(r2+d/2)); ln1 = log((r4-d/2)/(r3+d/2));
    Z3 = (eta/(2*pi))*ln3 * (1 - 0.5j*(Rs*(1/(r1+d/2)+1/(r2-d/2)))/(w*mu0*ln3));
    Z2 = (eta/(2*pi))*ln2 * (1 - 0.5j*(Rs*(1/(r2+d/2)+1/(r3-d/2)))/(w*mu0*ln2));
    Z1 = (eta/(2*pi))*ln1 * (1 - 0.5j*(Rs*(1/(r3+d/2)+1/(r4-d/2)))/(w*mu0*ln1));
    % 传播项
    a3 = 0.5*Rs/(eta*ln3)*(1/(r1+d/2)+1/(r2-d/2));
    a2 = 0.5*Rs/(eta*ln2)*(1/(r2+d/2)+1/(r3-d/2));
    a1 = 0.5*Rs/(eta*ln1)*(1/(r3+d/2)+1/(r4-d/2));
    b = w/c;
    t = @(a) (tanh(a*l)+1j*tan(b*l))/(1+1j*tan(b*l)*tanh(a*l));
    % 递归
    Zin1 = Z1*t(a1);
    Zin2 = Z2*(Zin1+Z2*t(a2))/(Z2+Zin1*t(a2));
    Zin = Z3*(Zin2+Z3*t(a3))/(Z3+Zin2*t(a3));
    ImY = imag(1/Zin);
end

