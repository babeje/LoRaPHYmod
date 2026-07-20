% plot_final_figures.m
% Построение четырёх финальных графиков для статьи.
% Требует двух файлов результатов:
%   r45 — прогон CR=4/5  (network_analysis_20260224_185118.mat)
%   r48 — прогон CR=4/8  (network_analysis_20260224_195225.mat)
%
% Сохраняет: fig3_per_vs_snr.fig, fig4_ber_vs_snr.fig,
%            fig5_e2e_per_vs_K.fig, fig6_e2e_thr_vs_K.fig

clear; clc;

% --- пути к файлам результатов ---
file_cr45 = fullfile('results', 'data', 'network_analysis_20260224_185118.mat');
file_cr48 = fullfile('results', 'data', 'network_analysis_20260224_195225.mat');

% --- папка для сохранения графиков ---
out_dir = fullfile('results', 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% --- загрузка данных ---
r45 = load(file_cr45);
r48 = load(file_cr48);

% --- общие параметры оформления ---
C_AWGN = [0.13, 0.47, 0.71];   % синий
C_TDL  = [0.17, 0.63, 0.17];   % зелёный
C_45   = [0.84, 0.15, 0.16];   % красный  — CR=4/5
C_48   = [0.12, 0.47, 0.71];   % синий    — CR=4/8
LW     = 1.6;
MS     = 6;
FS_ax  = 11;   % fontsize осей
FS_tt  = 12;   % fontsize заголовка
FS_leg = 9;

SNR_THRESH = -7.5;
FLOOR_VAL  = 1e-4;   % минимум для логарифмической оси

% вспомогательная функция: заменяем нули на FLOOR_VAL
clip = @(x) max(x, FLOOR_VAL);

% =========================================================
% Fig 3 — PER vs SNR (CR=4/5 и CR=4/8, AWGN и Rayleigh TDL)
% =========================================================
fig3 = figure('Name', 'PER vs SNR', 'NumberTitle', 'off', ...
              'Position', [100 100 640 440]);

semilogy(r45.snr_sweep, clip(r45.PER_awgn), ...
    '-o', 'Color', C_AWGN, 'LineWidth', LW, 'MarkerSize', MS, ...
    'DisplayName', 'AWGN, CR=4/5'); hold on;

semilogy(r48.snr_sweep, clip(r48.PER_awgn), ...
    '--^', 'Color', C_AWGN, 'LineWidth', LW, 'MarkerSize', MS, ...
    'DisplayName', 'AWGN, CR=4/8');

semilogy(r45.snr_sweep, clip(r45.PER_tdl), ...
    '-o', 'Color', C_TDL, 'LineWidth', LW, 'MarkerSize', MS, ...
    'DisplayName', 'Rayleigh TDL, CR=4/5');

semilogy(r48.snr_sweep, clip(r48.PER_tdl), ...
    '--^', 'Color', C_TDL, 'LineWidth', LW, 'MarkerSize', MS, ...
    'DisplayName', 'Rayleigh TDL, CR=4/8');

xline(SNR_THRESH, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, ...
    'Label', sprintf('SF7 threshold = %.1f dB', SNR_THRESH), ...
    'LabelOrientation', 'aligned', 'FontSize', 8);

xlabel('SNR, dB', 'FontSize', FS_ax);
ylabel('PER', 'FontSize', FS_ax);
title('PER vs SNR  |  SF=7, BW=125 kHz', 'FontSize', FS_tt);
legend('Location', 'southwest', 'FontSize', FS_leg);
xlim([r48.snr_sweep(1)-0.5, r48.snr_sweep(end)+0.5]);
ylim([5e-5, 2]);
grid on; grid minor;
set(gca, 'GridAlpha', 0.35, 'MinorGridAlpha', 0.15);

savefig(fig3, fullfile(out_dir, 'fig3_per_vs_snr.fig'));
fprintf('Saved: fig3_per_vs_snr.fig\n');

% =========================================================
% Fig 4 — BER vs SNR
% =========================================================
fig4 = figure('Name', 'BER vs SNR', 'NumberTitle', 'off', ...
              'Position', [120 120 640 440]);

semilogy(r45.snr_sweep, clip(r45.BER_awgn), ...
    '-o', 'Color', C_AWGN, 'LineWidth', LW, 'MarkerSize', MS, ...
    'DisplayName', 'AWGN, CR=4/5'); hold on;

semilogy(r48.snr_sweep, clip(r48.BER_awgn), ...
    '--^', 'Color', C_AWGN, 'LineWidth', LW, 'MarkerSize', MS, ...
    'DisplayName', 'AWGN, CR=4/8');

semilogy(r45.snr_sweep, clip(r45.BER_tdl), ...
    '-o', 'Color', C_TDL, 'LineWidth', LW, 'MarkerSize', MS, ...
    'DisplayName', 'Rayleigh TDL, CR=4/5');

semilogy(r48.snr_sweep, clip(r48.BER_tdl), ...
    '--^', 'Color', C_TDL, 'LineWidth', LW, 'MarkerSize', MS, ...
    'DisplayName', 'Rayleigh TDL, CR=4/8');

xline(SNR_THRESH, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, ...
    'Label', sprintf('SF7 threshold = %.1f dB', SNR_THRESH), ...
    'LabelOrientation', 'aligned', 'FontSize', 8);

xlabel('SNR, dB', 'FontSize', FS_ax);
ylabel('BER', 'FontSize', FS_ax);
title('BER vs SNR  |  SF=7, BW=125 kHz', 'FontSize', FS_tt);
legend('Location', 'southwest', 'FontSize', FS_leg);
xlim([r48.snr_sweep(1)-0.5, r48.snr_sweep(end)+0.5]);
ylim([5e-5, 1]);
grid on; grid minor;
set(gca, 'GridAlpha', 0.35, 'MinorGridAlpha', 0.15);

savefig(fig4, fullfile(out_dir, 'fig4_ber_vs_snr.fig'));
fprintf('Saved: fig4_ber_vs_snr.fig\n');

% =========================================================
% Fig 5 — E2E PER vs K
% =========================================================
K45 = r45.hop_range;
K48 = r48.hop_range;
snr_typical = r48.snr_typical;   % worst-case SNR_hop

fig5 = figure('Name', 'E2E PER vs K', 'NumberTitle', 'off', ...
              'Position', [140 140 640 420]);

plot(K45, r45.PER_e2e_hops, '-o', 'Color', C_45, ...
    'LineWidth', LW, 'MarkerSize', MS, 'DisplayName', 'CR=4/5'); hold on;

plot(K48, r48.PER_e2e_hops, '--s', 'Color', C_48, ...
    'LineWidth', LW, 'MarkerSize', MS, 'DisplayName', 'CR=4/8');

xline(7, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, ...
    'Label', 'Route (7 hops)', 'LabelOrientation', 'aligned', 'FontSize', 8);

xlabel('Number of hops K', 'FontSize', FS_ax);
ylabel('PER_{e2e}', 'FontSize', FS_ax);
title(sprintf('End-to-end PER vs Hops K  |  worst-case SNR_{hop} = %.1f dB', ...
    snr_typical), 'FontSize', FS_tt);
legend('Location', 'northwest', 'FontSize', FS_leg);
xlim([0.5, max(K45(end), K48(end))+0.5]);
ylim([0, 1]);
grid on;
set(gca, 'GridAlpha', 0.35);

savefig(fig5, fullfile(out_dir, 'fig5_e2e_per_vs_K.fig'));
fprintf('Saved: fig5_e2e_per_vs_K.fig\n');

% =========================================================
% Fig 6 — E2E Throughput vs K
% =========================================================
fig6 = figure('Name', 'E2E Throughput vs K', 'NumberTitle', 'off', ...
              'Position', [160 160 640 420]);

plot(K45, r45.Thr_e2e_hops, '-o', 'Color', C_45, ...
    'LineWidth', LW, 'MarkerSize', MS, 'DisplayName', 'CR=4/5'); hold on;

plot(K48, r48.Thr_e2e_hops, '--s', 'Color', C_48, ...
    'LineWidth', LW, 'MarkerSize', MS, 'DisplayName', 'CR=4/8');

xline(7, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, ...
    'Label', 'Route (7 hops)', 'LabelOrientation', 'aligned', 'FontSize', 8);

xlabel('Number of hops K', 'FontSize', FS_ax);
ylabel('Throughput_{e2e}, bit/s', 'FontSize', FS_ax);
title(sprintf('End-to-end Throughput vs Hops K  |  worst-case SNR_{hop} = %.1f dB', ...
    snr_typical), 'FontSize', FS_tt);
legend('Location', 'northeast', 'FontSize', FS_leg);
xlim([0.5, max(K45(end), K48(end))+0.5]);
ylim([0, max(r45.Thr_e2e_hops(1), r48.Thr_e2e_hops(1)) * 1.1]);
grid on;
set(gca, 'GridAlpha', 0.35);

savefig(fig6, fullfile(out_dir, 'fig6_e2e_thr_vs_K.fig'));
fprintf('Saved: fig6_e2e_thr_vs_K.fig\n');

fprintf('\nВсе графики сохранены в: %s\n', out_dir);
