clear; clc; close all;

projectRoot = fileparts(mfilename('fullpath'));
projectRoot = fileparts(projectRoot);
addpath(genpath(projectRoot));

%% ============================================================
%  PHY-параметры
%% ============================================================
rf_freq = 868e6;
sf      = 7;
bw      = 125e3;
fs      = 1e6;

CR          = 1;        % code rate 4/5 — соответствует snrThreshold=-7.5 дБ
HasHeader   = true;
UseCRC      = true;
PreambleLen = 8;
FastMode    = true;

Npkts          = 100;
payloadLenBits = 128;

snr_list = -12:1:5;

% Порог чувствительности SF7 по спецификации SX1262
snrThreshold = -7.5;

%% ============================================================
%  Сценарии: sweep по скорости (Doppler) с TDL-каналом
%% ============================================================
cases(1).name  = "v=0 м/с (статика)";
cases(1).v_mps = 0;

cases(2).name  = "v=5 м/с (пешеход)";
cases(2).v_mps = 5;

cases(3).name  = "v=15 м/с (UGV)";
cases(3).v_mps = 15;

cases(4).name  = "v=30 м/с (БПЛА)";
cases(4).v_mps = 30;

% TDL-профиль (4 луча)
tdlDelays = [0, 0.5e-6, 2.0e-6, 4.0e-6];
tdlGains  = [0, -4,     -10,    -18];

PER = zeros(numel(cases), numel(snr_list));
BER = zeros(numel(cases), numel(snr_list));

%% ============================================================
%  Основной цикл
%% ============================================================
for ci = 1:numel(cases)
    v_mps  = cases(ci).v_mps;
    fd_Hz  = (v_mps / 3e8) * rf_freq;   % f_D = v/c * f_c

    fprintf('\n=== %s (f_D = %.1f Гц) ===\n', cases(ci).name, fd_Hz);

    for si = 1:numel(snr_list)
        snr_dB = snr_list(si);

        modem = LoRaModem(rf_freq, sf, bw, fs, ...
            'CR', CR, 'HasHeader', HasHeader, 'UseCRC', UseCRC, ...
            'PreambleLen', PreambleLen, 'FastMode', FastMode);

        % Seed=[] — статистически независимые реализации для каждой
        % точки SNR. Усреднение по Npkts пакетам обеспечивает
        % сходимость к истинному E{PER}(SNR).
        channel = RayleighTDLChannel(fs, snr_dB, 0, ...
            'SF',         sf, ...
            'BW',         bw, ...
            'DopplerHz',  fd_Hz, ...
            'PathDelays', tdlDelays, ...
            'PathGains',  tdlGains, ...
            'Seed',       []);

        sim = LoRaSimulator(modem, channel);
        [ber, per] = sim.run(Npkts, payloadLenBits);

        BER(ci, si) = ber;
        PER(ci, si) = per;

        fprintf('  SNR=%5.1f дБ | BER=%9.3e | PER=%6.4f\n', snr_dB, ber, per);
    end
end

%% ============================================================
%  Графики
%% ============================================================
colors = [0.00 0.45 0.70;
          0.47 0.67 0.19;
          0.85 0.33 0.10;
          0.63 0.08 0.18];

% --- PER vs SNR ---
figure('Name', 'PER vs SNR', 'Color', 'w', 'Position', [50 50 720 500]);
hold on;
for ci = 1:numel(cases)
    semilogy(snr_list, max(PER(ci,:), 1e-4), '-o', ...
        'Color', colors(ci,:), 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerFaceColor', colors(ci,:));
end
xline(snrThreshold, '--k', 'LineWidth', 1.3, ...
    'Label', sprintf('SF%d threshold = %.1f dB', sf, snrThreshold), ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 10);
xlabel('SNR, дБ', 'FontSize', 13);
ylabel('PER', 'FontSize', 13);
title(sprintf('PER vs SNR | SF=%d, BW=%d кГц, CR=4/%d | Rayleigh TDL + Doppler', ...
    sf, bw/1e3, CR+4), 'FontSize', 13);
legend({cases.name}, 'Location', 'southwest', 'FontSize', 11);
grid on;
ylim([1e-3 1]);
xlim([snr_list(1) snr_list(end)]);
hold off;

% --- BER vs SNR ---
figure('Name', 'BER vs SNR', 'Color', 'w', 'Position', [70 70 720 500]);
hold on;
for ci = 1:numel(cases)
    semilogy(snr_list, max(BER(ci,:), 1e-4), '-o', ...
        'Color', colors(ci,:), 'LineWidth', 2, 'MarkerSize', 6, ...
        'MarkerFaceColor', colors(ci,:));
end
xline(snrThreshold, '--k', 'LineWidth', 1.3, ...
    'Label', sprintf('SF%d threshold = %.1f dB', sf, snrThreshold), ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 10);
xlabel('SNR, дБ', 'FontSize', 13);
ylabel('BER', 'FontSize', 13);
title(sprintf('BER vs SNR | SF=%d, BW=%d кГц, CR=4/%d | Rayleigh TDL + Doppler', ...
    sf, bw/1e3, CR+4), 'FontSize', 13);
legend({cases.name}, 'Location', 'southwest', 'FontSize', 11);
grid on;
ylim([1e-3 1]);
xlim([snr_list(1) snr_list(end)]);
hold off;
% % --- PER vs SNR ---
% figure('Name', 'PER vs SNR', 'Color', 'w', 'Position', [50 50 720 500]);
% hold on;
% for ci = 1:numel(cases)
%     semilogy(snr_list, max(PER(ci,:), 1e-4), '-o', ...
%         'Color', colors(ci,:), 'LineWidth', 2, 'MarkerSize', 6, ...
%         'MarkerFaceColor', colors(ci,:));
% end
% xline(snrThreshold, '--k', 'LineWidth', 1.3, ...
%     'Label', sprintf('SF%d threshold = %.1f dB', sf, snrThreshold), ...
%     'LabelHorizontalAlignment', 'left', 'FontSize', 10);
% xlabel('SNR, дБ', 'FontSize', 13);
% ylabel('PER', 'FontSize', 13);
% title(sprintf('PER vs SNR | SF=%d, BW=%d кГц, CR=4/%d | Rayleigh TDL + Doppler', ...
%     sf, bw/1e3, CR+4), 'FontSize', 13);
% legend({cases.name}, 'Location', 'southwest', 'FontSize', 11);
% grid on; ylim([1e-4 1]); xlim([snr_list(1) snr_list(end)]);
% hold off;
% 
% % --- BER vs SNR ---
% figure('Name', 'BER vs SNR', 'Color', 'w', 'Position', [70 70 720 500]);
% hold on;
% for ci = 1:numel(cases)
%     semilogy(snr_list, max(BER(ci,:), 1e-5), '-o', ...
%         'Color', colors(ci,:), 'LineWidth', 2, 'MarkerSize', 6, ...
%         'MarkerFaceColor', colors(ci,:));
% end
% xline(snrThreshold, '--k', 'LineWidth', 1.3, ...
%     'Label', sprintf('SF%d threshold = %.1f dB', sf, snrThreshold), ...
%     'LabelHorizontalAlignment', 'left', 'FontSize', 10);
% xlabel('SNR, дБ', 'FontSize', 13);
% ylabel('BER', 'FontSize', 13);
% title(sprintf('BER vs SNR | SF=%d, BW=%d кГц, CR=4/%d | Rayleigh TDL + Doppler', ...
%     sf, bw/1e3, CR+4), 'FontSize', 13);
% legend({cases.name}, 'Location', 'southwest', 'FontSize', 11);
% grid on; ylim([1e-5 1]); xlim([snr_list(1) snr_list(end)]);
% hold off;

%% ============================================================
%  Сохранение
%% ============================================================
outDir = fullfile(projectRoot, 'results', 'data');
if ~exist(outDir, 'dir'), mkdir(outDir); end

ts = datestr(now, 'yyyymmdd_HHMMSS');
save(fullfile(outDir, "snr_sweep_" + ts + ".mat"), ...
    'snr_list', 'cases', 'BER', 'PER', ...
    'rf_freq', 'sf', 'bw', 'fs', 'Npkts', 'payloadLenBits', ...
    'CR', 'HasHeader', 'UseCRC', 'PreambleLen', 'FastMode', ...
    'tdlDelays', 'tdlGains', 'snrThreshold');

fprintf('\nРезультаты сохранены: results/data/snr_sweep_%s.mat\n', ts);