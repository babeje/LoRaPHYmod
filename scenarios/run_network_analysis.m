%% run_network_analysis.m
% Сценарий: оценка энергоэффективности и пропускной способности
% децентрализованной LoRa-сети с учётом динамически меняющихся условий.
%
% Структура сценария:
%   Блок 1  — Параметры сети и топология
%   Блок 2  — PHY-параметры LoRa
%   Блок 3  — SNR-матрица по ITU-R P.1411
%   Блок 4  — Маршрутизация и SNR по хопам
%   Блок 5  — PHY-симуляция: BER/PER по каждому хопу маршрута
%   Блок 6  — End-to-end метрики: PER_e2e, Throughput_e2e
%   Блок 7  — Кривые BER/PER vs SNR (один хоп, три канала)
%   Блок 8  — PER_e2e и Throughput vs число хопов (аналитика)
%   Блок 9  — Влияние скорости движения (DopplerChannel)
%   Блок 10 — Графики

clear; clc; close all;

projectRoot = fileparts(mfilename('fullpath'));
projectRoot = fileparts(projectRoot);
addpath(genpath(projectRoot));

fprintf('=== Анализ LoRa mesh-сети: PHY-уровень ===\n\n');

%% ============================================================
%  БЛОК 1: ПАРАМЕТРЫ СЕТИ И ТОПОЛОГИЯ
%% ============================================================
% Параметры топологии подобраны под многохоповый маршрут:
%
%   areaSize = 1200 м — расчётное значение из модели распространения.
%                       При txPower=10 дБм и NLOS с экспонентой пути n=3.8
%                       максимальный радиус одного хопа составляет ~520 м.
%                       Диагональ площадки = 1200*sqrt(2) ≈ 1700 м ≈ 3*d_hop.
%                       Это обеспечивает 3-4 хопа на маршруте src→dst.
%
%   Важно: при старой NLOS-формуле (L_excess = 20 + 30*log10(d_km))
%          потери при d < 1 км были занижены (отрицательный L_excess),
%          из-за чего все узлы видели друг друга напрямую и маршрут
%          вырождался в 1 хоп независимо от areaSize. Проблема устранена
%          в ituP1411_corrected: NLOS теперь использует степенной закон
%          L = L_ref(1 м) + 10 * n_path * log10(d), n_path = 3.8.
%
%   nodeHeights 5..25 м — реалистичный монтаж на зданиях/столбах.
%
%   numNodes = 30     — достаточно для связной топологии без избыточной
%                       плотности, коллапсирующей маршрут в 1 хоп.

numNodes = 30;
areaSize = 1200;  % м — подобрано из условия: d_hop_max ≈ 520 м → 3-4 хопа

rng(42);
nodePositions = rand(numNodes, 2) * areaSize;

rng(42);
nodeHeights = 5 + rand(numNodes, 1) * 20;   % 5..25 м

X = nodePositions(:, 1);
Y = nodePositions(:, 2);
Z = nodeHeights;

distances  = zeros(numNodes);
elevAngles = zeros(numNodes);
for i = 1:numNodes
    for j = 1:numNodes
        dx = X(i) - X(j);
        dy = Y(i) - Y(j);
        dz = Z(i) - Z(j);
        distances(i,j)  = sqrt(dx^2 + dy^2 + dz^2);
        elevAngles(i,j) = atan2d(dz, sqrt(dx^2 + dy^2));
    end
end

%% ============================================================
%  БЛОК 2: PHY-ПАРАМЕТРЫ LoRa
%% ============================================================
% Параметры соответствуют EBYTE E22-868T22U (SX1262, EU868).
%
%   txPower_dBm = 10  — намеренно снижено относительно максимума (22 дБм),
%                       чтобы радиус хопа составил ~80-120 м в NLOS.
%                       При areaSize=300 м это обеспечивает 3-5 хопов
%                       на маршруте между наиболее удалёнными узлами.

fc_MHz       = 868;
fc_Hz        = fc_MHz * 1e6;
txPower_dBm  = 10;
noiseFigure  = 6;           % дБ, типовой SX1262
hTx_m        = 1.5;         % высота антенны над точкой монтажа, м
hRx_m        = 1.5;
snrThreshold = -7.5;        % дБ, порог SF7

sf          = 7;
bw          = 125e3;
fs          = 1e6;
CR          = 1;
payloadBits = 128;          % 16 байт — типичная телеметрия UAV/UGV
Npkts       = 1000;
Npkts_sweep = 200;

k_B            = 1.38e-23;
T_K            = 290;
noisePower_dBm = 10*log10(k_B * T_K * bw * 1000) + noiseFigure;

T_sym_ms  = (2^sf / bw) * 1e3;
Nsym_pkt  = 8 + ceil((payloadBits/8) * 2 / (sf-2)) * (CR+4);
T_pkt_ms  = (8 + 4.25 + Nsym_pkt) * T_sym_ms;
T_pkt_s   = T_pkt_ms / 1e3;

% TDL-профиль: 4 луча, tau_rms ≈ 1.75 мкс (ITU-R P.1411-10, 300 м)
tdlDelays = [0, 0.5e-6, 2.0e-6, 5.0e-6];
tdlGains  = [0, -2,     -5,     -8];

snr_sweep = -12:1:5;

fprintf('Параметры PHY:\n');
fprintf('  SF=%d, BW=%.0f кГц, CR=4/%d, payload=%d бит\n', sf, bw/1e3, CR+4, payloadBits);
fprintf('  T_sym=%.3f мс,  T_pkt≈%.1f мс\n', T_sym_ms, T_pkt_ms);
fprintf('  P_noise=%.1f дБм,  SNR_thresh=%.1f дБ\n\n', noisePower_dBm, snrThreshold);

%% ============================================================
%  БЛОК 3: SNR-МАТРИЦА (ITU-R P.1411)
%% ============================================================

pathLoss  = zeros(numNodes);
snrMatrix = zeros(numNodes);

for i = 1:numNodes
    for j = 1:numNodes
        if i ~= j
            hTx_eff = Z(i) + hTx_m;
            hRx_eff = Z(j) + hRx_m;
            L = ituP1411_corrected(fc_MHz, distances(i,j), ...
                hTx_eff, hRx_eff, elevAngles(i,j));
            pathLoss(i,j)  = L;
            snrMatrix(i,j) = txPower_dBm - L - noisePower_dBm;
        else
            snrMatrix(i,j) = NaN;
        end
    end
end

connectivity   = snrMatrix > snrThreshold;
connectedNodes = sum(any(connectivity, 2));
snrConn        = snrMatrix(connectivity);

fprintf('Топология:\n');
fprintf('  Связанных узлов: %d из %d\n', connectedNodes, numNodes);
fprintf('  SNR: min=%.1f  median=%.1f  max=%.1f дБ\n\n', ...
    min(snrConn), median(snrConn), max(snrConn));

%% ============================================================
%  БЛОК 4: МАРШРУТИЗАЦИЯ И SNR ПО ХОПАМ
%% ============================================================

[src, dst] = findMostDistantNodes(nodePositions);
route      = findRoute(connectivity, src, dst);

% Если маршрут слишком короткий — ищем пару с максимальным числом хопов
if numel(route) < 4
    bestLen = numel(route);
    for s = 1:numNodes
        for d = 1:numNodes
            if s ~= d
                r = findRoute(connectivity, s, d);
                if numel(r) > bestLen
                    bestLen = numel(r);
                    src = s; dst = d; route = r;
                end
            end
        end
    end
    fprintf('Найден более длинный маршрут: %d→%d (%d хопов)\n\n', ...
        src, dst, numel(route)-1);
end

if isempty(route)
    error('Маршрут не найден. Увеличьте areaSize или numNodes.');
end

nHops     = length(route) - 1;
routeSNR  = zeros(nHops, 1);
routeDist = zeros(nHops, 1);

fprintf('Маршрут: узел %d → %d  (%d хопов)\n', src, dst, nHops);
for h = 1:nHops
    n1 = route(h);  n2 = route(h+1);
    routeSNR(h)  = snrMatrix(n1, n2);
    routeDist(h) = distances(n1, n2);
    fprintf('  Хоп %d: %2d→%2d  d=%5.0f м  SNR=%5.1f дБ\n', ...
        h, n1, n2, routeDist(h), routeSNR(h));
end
fprintf('\n');

%% ============================================================
%  БЛОК 5: PHY-СИМУЛЯЦИЯ ПО ХОПАМ МАРШРУТА
%% ============================================================

fprintf('PHY-симуляция по хопам (Npkts=%d)...\n', Npkts);

hopBER     = nan(nHops, 1);
hopPER     = nan(nHops, 1);
hopThr_bps = nan(nHops, 1);

for h = 1:nHops
    modem = LoRaModem(fc_Hz, sf, bw, fs, ...
        'CR', CR, 'HasHeader', true, 'UseCRC', true, ...
        'PreambleLen', 8, 'FastMode', true);

    channel = RayleighTDLChannel(fs, routeSNR(h), 0, ...
        'PathDelays', tdlDelays, ...
        'PathGains',  tdlGains, ...
        'Seed',       42 + h);

    [hopBER(h), hopPER(h), ~] = LoRaSimulator(modem, channel).run(Npkts, payloadBits);
    hopThr_bps(h) = (1 - hopPER(h)) * payloadBits / T_pkt_s;

    fprintf('  Хоп %d: SNR=%5.1f дБ | BER=%.2e | PER=%.4f | Thr=%6.0f бит/с\n', ...
        h, routeSNR(h), hopBER(h), hopPER(h), hopThr_bps(h));
end

%% ============================================================
%  БЛОК 6: END-TO-END МЕТРИКИ
%% ============================================================

PER_e2e         = 1 - prod(1 - hopPER);
Thr_e2e         = min(hopThr_bps);
delay_no_arq_ms = nHops * T_pkt_ms;
delay_arq_ms    = sum(T_pkt_ms ./ max(1 - hopPER, 1e-9));

fprintf('\nEnd-to-end:\n');
fprintf('  PER_e2e          = %.4f  (%.1f%%)\n', PER_e2e, PER_e2e*100);
fprintf('  Throughput_e2e   = %.0f бит/с\n',    Thr_e2e);
fprintf('  Задержка (нет ARQ) = %.1f мс\n',     delay_no_arq_ms);
fprintf('  Задержка (ARQ)     = %.1f мс\n\n',   delay_arq_ms);

%% ============================================================
%  БЛОК 7: BER/PER vs SNR — ОДИН ХОП, ТРИ КАНАЛА
%% ============================================================

fprintf('Sweep BER/PER vs SNR (%d точек, Npkts=%d)...\n', numel(snr_sweep), Npkts_sweep);

BER_awgn     = nan(size(snr_sweep));  BER_rayleigh = nan(size(snr_sweep));  BER_tdl = nan(size(snr_sweep));
PER_awgn     = nan(size(snr_sweep));  PER_rayleigh = nan(size(snr_sweep));  PER_tdl = nan(size(snr_sweep));

for si = 1:numel(snr_sweep)
    snr_i = snr_sweep(si);

    modem = LoRaModem(fc_Hz, sf, bw, fs, ...
        'CR', CR, 'HasHeader', true, 'UseCRC', true, ...
        'PreambleLen', 8, 'FastMode', true);

    [BER_awgn(si), PER_awgn(si), ~] = ...
        LoRaSimulator(modem, AwgnChannel(fs, snr_i)).run(Npkts_sweep, payloadBits);

    ch_ray = RayleighTDLChannel(fs, snr_i, 0, 'PathDelays', 0, 'PathGains', 0, 'Seed', 42);
    [BER_rayleigh(si), PER_rayleigh(si), ~] = ...
        LoRaSimulator(modem, ch_ray).run(Npkts_sweep, payloadBits);

    ch_tdl = RayleighTDLChannel(fs, snr_i, 0, ...
        'PathDelays', tdlDelays, 'PathGains', tdlGains, 'Seed', 42);
    [BER_tdl(si), PER_tdl(si), ~] = ...
        LoRaSimulator(modem, ch_tdl).run(Npkts_sweep, payloadBits);

    fprintf('  SNR=%5.1f дБ | PER: AWGN=%.2f  Ray=%.2f  TDL=%.2f\n', ...
        snr_i, PER_awgn(si), PER_rayleigh(si), PER_tdl(si));
end

%% ============================================================
%  БЛОК 8: PER_E2E И THROUGHPUT VS ЧИСЛО ХОПОВ
%% ============================================================

fprintf('\nАналитика E2E vs K...\n');

snr_typical = median(routeSNR);
ch_h = RayleighTDLChannel(fs, snr_typical, 0, ...
    'PathDelays', tdlDelays, 'PathGains', tdlGains, 'Seed', 42);
modem_h = LoRaModem(fc_Hz, sf, bw, fs, ...
    'CR', CR, 'HasHeader', true, 'UseCRC', true, 'PreambleLen', 8, 'FastMode', true);
[~, PER_single, ~] = LoRaSimulator(modem_h, ch_h).run(Npkts, payloadBits);

hop_range    = 1:12;
PER_e2e_hops = 1 - (1 - PER_single).^hop_range;
Thr_e2e_hops = (1 - PER_single) * payloadBits ./ (hop_range * T_pkt_s);

fprintf('  SNR_typical=%.1f дБ,  PER_hop=%.4f\n', snr_typical, PER_single);

%% ============================================================
%  БЛОК 9: ВЛИЯНИЕ СКОРОСТИ ДВИЖЕНИЯ
%% ============================================================

fprintf('\nSweep по скоростям...\n');

velocities  = [0, 5, 15, 30];
PER_dynamic = nan(numel(velocities), numel(snr_sweep));

for vi = 1:numel(velocities)
    v  = velocities(vi);
    fd = (v / 3e8) * fc_Hz;
    for si = 1:numel(snr_sweep)
        modem_d = LoRaModem(fc_Hz, sf, bw, fs, ...
            'CR', CR, 'HasHeader', true, 'UseCRC', true, ...
            'PreambleLen', 8, 'FastMode', true);
        ch_d = DopplerChannel(fs, snr_sweep(si), 0, fc_Hz, v, 0, 0);
        [~, PER_dynamic(vi, si), ~] = ...
            LoRaSimulator(modem_d, ch_d).run(Npkts_sweep, payloadBits);
    end
    idx0 = find(snr_sweep == 0, 1);
    if ~isempty(idx0)
        fprintf('  v=%2d м/с (fd=%.0f Гц): PER@0дБ=%.3f\n', v, fd, PER_dynamic(vi, idx0));
    end
end

%% ============================================================
%  БЛОК 10: ГРАФИКИ
%% ============================================================

colors = [0.00 0.45 0.70;
          0.85 0.33 0.10;
          0.47 0.67 0.19;
          0.63 0.08 0.18];

% --- График 1: PER vs SNR ---
figure('Name','PER vs SNR','Color','w','Position',[50 50 680 480]);
semilogy(snr_sweep, max(PER_awgn,     1e-4), '-o', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 6); hold on;
semilogy(snr_sweep, max(PER_rayleigh, 1e-4), '-s', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 6);
semilogy(snr_sweep, max(PER_tdl,      1e-4), '-^', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 6);
xline(snrThreshold, '--k', 'LineWidth', 1.2, 'Label', sprintf('Порог SF%d', sf));
xlabel('SNR, дБ', 'FontSize', 13); ylabel('PER', 'FontSize', 13);
title(sprintf('PER vs SNR | SF=%d, BW=%d кГц, CR=4/%d', sf, bw/1e3, CR+4), 'FontSize', 13);
legend('AWGN','Rayleigh (плоский)','Rayleigh TDL (городской)','Location','southwest');
grid on; ylim([1e-3 1]); xlim([snr_sweep(1) snr_sweep(end)]);

% --- График 2: BER vs SNR ---
figure('Name','BER vs SNR','Color','w','Position',[70 70 680 480]);
semilogy(snr_sweep, max(BER_awgn,     1e-5), '-o', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 6); hold on;
semilogy(snr_sweep, max(BER_rayleigh, 1e-5), '-s', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 6);
semilogy(snr_sweep, max(BER_tdl,      1e-5), '-^', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 6);
xlabel('SNR, дБ', 'FontSize', 13); ylabel('BER', 'FontSize', 13);
title(sprintf('BER vs SNR | SF=%d, BW=%d кГц, CR=4/%d', sf, bw/1e3, CR+4), 'FontSize', 13);
legend('AWGN','Rayleigh (плоский)','Rayleigh TDL (городской)','Location','southwest');
grid on; ylim([1e-4 1]); xlim([snr_sweep(1) snr_sweep(end)]);

% --- График 3: SNR и PER по хопам маршрута ---
if nHops > 1
    tickLbls = arrayfun(@(h) sprintf('%d→%d', route(h), route(h+1)), 1:nHops, 'UniformOutput', false);
    figure('Name','Хопы маршрута','Color','w','Position',[90 90 760 430]);
    yyaxis left;
    bar(1:nHops, routeSNR, 0.55, 'FaceColor', colors(1,:), 'FaceAlpha', 0.75);
    ylabel('SNR хопа, дБ', 'FontSize', 13);
    ylim([min(routeSNR)-5, max(routeSNR)+5]);
    yyaxis right;
    plot(1:nHops, hopPER, '-s', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 9);
    ylabel('PER хопа', 'FontSize', 13); ylim([0 1]);
    set(gca, 'XTick', 1:nHops, 'XTickLabel', tickLbls, 'XTickLabelRotation', 35);
    xlabel('Хоп', 'FontSize', 13);
    title(sprintf('Маршрут %d→%d | %d хопов | PER_{e2e}=%.3f', src, dst, nHops, PER_e2e), 'FontSize', 13);
    legend('SNR хопа','PER хопа','Location','best'); grid on;
end

% --- График 4: E2E метрики vs число хопов ---
figure('Name','E2E vs хопы','Color','w','Position',[110 110 720 430]);
yyaxis left;
plot(hop_range, PER_e2e_hops, '-o', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8);
ylabel('PER_{e2e}', 'FontSize', 13); ylim([0 1]);
yyaxis right;
plot(hop_range, Thr_e2e_hops, '-s', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Throughput_{e2e}, бит/с', 'FontSize', 13);
xline(nHops, '--k', 'LineWidth', 1.2, 'Label', sprintf('Маршрут (%d хопов)', nHops));
xlabel('Число хопов K', 'FontSize', 13);
title(sprintf('E2E метрики vs K | SNR_{хоп}=%.1f дБ', snr_typical), 'FontSize', 13);
legend('PER_{e2e}','Throughput_{e2e}','Location','east'); grid on; xlim([1 hop_range(end)]);

% --- График 5: PER vs SNR при разных скоростях ---
figure('Name','PER vs SNR — скорость','Color','w','Position',[130 130 680 480]);
vel_labels = cell(numel(velocities), 1);
for vi = 1:numel(velocities)
    fd_v = (velocities(vi) / 3e8) * fc_Hz;
    semilogy(snr_sweep, max(PER_dynamic(vi,:), 1e-4), '-', 'Color', colors(vi,:), 'LineWidth', 2); hold on;
    vel_labels{vi} = sprintf('v = %d м/с  (f_D = %.0f Гц)', velocities(vi), fd_v);
end
xlabel('SNR, дБ', 'FontSize', 13); ylabel('PER', 'FontSize', 13);
title(sprintf('PER vs SNR | SF=%d, влияние скорости', sf), 'FontSize', 13);
legend(vel_labels, 'Location', 'southwest');
grid on; ylim([1e-3 1]); xlim([snr_sweep(1) snr_sweep(end)]);

% --- График 6: 3D-сеть ---
figure('Name','3D Network','Color','w','Position',[150 150 800 580]);
hold on;
for i = 1:numNodes
    for j = i+1:numNodes
        if connectivity(i,j)
            plot3([X(i) X(j)],[Y(i) Y(j)],[Z(i) Z(j)], ...
                'Color',[0.75 0.75 0.75],'LineWidth',0.5,'HandleVisibility','off');
        end
    end
end
for h = 1:nHops
    c = [hopPER(h), 1-hopPER(h), 0];
    plot3(X(route(h:h+1)), Y(route(h:h+1)), Z(route(h:h+1)), '-', 'Color', c, 'LineWidth', 4.5);
end
scatter3(X, Y, Z, 70, Z, 'filled');
colormap(jet); colorbar;
xlabel('X, м'); ylabel('Y, м'); zlabel('Высота, м');
title(sprintf('3D-сеть | маршрут %d→%d | цвет хопа = PER', src, dst), 'FontSize', 13);
grid on; box on; view(45, 30); rotate3d on; hold off;

%% ============================================================
%  СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
%% ============================================================

outDir = fullfile(projectRoot, 'results', 'data');
if ~exist(outDir, 'dir'), mkdir(outDir); end

ts    = datestr(now, 'yyyymmdd_HHMMSS');
fname = fullfile(outDir, ['network_analysis_' ts '.mat']);
save(fname, ...
    'snr_sweep', ...
    'BER_awgn','BER_rayleigh','BER_tdl', ...
    'PER_awgn','PER_rayleigh','PER_tdl', ...
    'routeSNR','routeDist','hopBER','hopPER','hopThr_bps', ...
    'PER_e2e','Thr_e2e','delay_no_arq_ms','delay_arq_ms', ...
    'PER_e2e_hops','Thr_e2e_hops','hop_range', ...
    'PER_dynamic','velocities', ...
    'route','src','dst','nHops', ...
    'sf','bw','payloadBits','Npkts','T_pkt_ms', ...
    'tdlDelays','tdlGains','snr_typical','PER_single');

fprintf('\nРезультаты сохранены: %s\n', fname);
fprintf('=== Симуляция завершена ===\n');
