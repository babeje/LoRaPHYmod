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
areaSize = 800;   % м — при txPower=20 дБм и d_hop~450 м → 3-4 хопа, SNR>5 дБ

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
txPower_dBm  = 14;          % 20 дБм → d_hop ~450 м при SNR>5 дБ (NLOS, n=3.8)
noiseFigure  = 6;           % дБ, типовой SX1262
hTx_m        = 1.5;         % высота антенны над точкой монтажа, м
hRx_m        = 1.5;
snrThreshold = -7.5;        % дБ, порог SF7

sf          = 7;
bw          = 125e3;
fs          = 1e6;
CR          = 4;
payloadBits = 128;          % 16 байт — типичная телеметрия UAV/UGV
Npkts       = 1000;
Npkts_sweep = 1000;

k_B            = 1.38e-23;
T_K            = 290;
noisePower_dBm = 10*log10(k_B * T_K * bw * 1000) + noiseFigure;

T_sym_ms  = (2^sf / bw) * 1e3;
Nsym_pkt  = 8 + ceil((payloadBits/8) * 2 / (sf-2)) * (CR+4);
T_pkt_ms  = (8 + 4.25 + Nsym_pkt) * T_sym_ms;
T_pkt_s   = T_pkt_ms / 1e3;

% TDL-профиль: 4 луча, tau_rms ≈ 1.75 мкс (ITU-R P.1411-10, 300 м)
PathDelays = [0,  0.5e-6,  2.0e-6,  4.0e-6];
PathGains  = [0,  -4,      -10,     -18];

snr_sweep = -12:1:10;     % расширен до +10 дБ — TDL успевает достичь PER≈0

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
        'PreambleLen', 8, 'FastMode', false);

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

    % FastMode=false: полное декодирование — точные кривые BER/PER
    modem = LoRaModem(fc_Hz, sf, bw, fs, ...
        'CR', CR, 'HasHeader', true, 'UseCRC', true, ...
        'PreambleLen', 8, 'FastMode', false);

    [BER_awgn(si), PER_awgn(si), ~] = ...
        LoRaSimulator(modem, AwgnChannel(fs, snr_i)).run(Npkts_sweep, payloadBits);

    ch_ray = RayleighTDLChannel(fs, snr_i, 0, 'PathDelays', 0, 'PathGains', 0, 'Seed', 42);
    [BER_rayleigh(si), PER_rayleigh(si), ~] = ...
        LoRaSimulator(modem, ch_ray).run(Npkts_sweep, payloadBits);

    ch_tdl = RayleighTDLChannel(fs, snr_i, 0, ...
    'PathDelays', tdlDelays, 'PathGains', tdlGains, ...
    'Seed', []);   % [] = случайный, усредняем по реализациям
    % ch_tdl = RayleighTDLChannel(fs, snr_i, 0, ...
    %     'PathDelays', tdlDelays, 'PathGains', tdlGains, 'Seed', 42);
    % [BER_tdl(si), PER_tdl(si), ~] = ...
    %     LoRaSimulator(modem, ch_tdl).run(Npkts_sweep, payloadBits);

    fprintf('  SNR=%5.1f дБ | PER: AWGN=%.2f  Ray=%.2f  TDL=%.2f\n', ...
        snr_i, PER_awgn(si), PER_rayleigh(si), PER_tdl(si));
end

%% ============================================================
%  БЛОК 8: PER_E2E И THROUGHPUT VS ЧИСЛО ХОПОВ
%% ============================================================

fprintf('\nАналитика E2E vs K...\n');

% Типовой SNR хопа — берём минимальный SNR маршрута (worst-case hop)
% Это честнее median: если хоть один хоп слабый — он определяет PER_e2e
snr_typical  = min(routeSNR);
snr_typical_median = median(routeSNR);
ch_h = RayleighTDLChannel(fs, snr_typical, 0, ...
    'PathDelays', tdlDelays, 'PathGains', tdlGains, 'Seed', 42);
% FastMode=false для достоверной оценки PER_hop на типовом SNR
modem_h = LoRaModem(fc_Hz, sf, bw, fs, ...
    'CR', CR, 'HasHeader', true, 'UseCRC', true, 'PreambleLen', 8, 'FastMode', false);
[~, PER_single, ~] = LoRaSimulator(modem_h, ch_h).run(Npkts, payloadBits);

hop_range    = 1:12;
% PER_e2e(K) = 1 - (1-PER_hop)^K  — растёт с K: чем больше хопов,
%              тем выше вероятность потери пакета хотя бы на одном из них.
PER_e2e_hops = 1 - (1 - PER_single).^hop_range;

% Throughput_e2e(K): полезная скорость на выходе маршрута.
% T_total = K * T_pkt — суммарное время передачи через K хопов.
% Доля успешных пакетов = (1-PER_hop)^K (все хопы должны пройти).
% Throughput = успешные биты / суммарное время:
%   Thr = (1-PER_hop)^K * payloadBits / (K * T_pkt_s)
% Числитель убывает как (1-PER_hop)^K, знаменатель растёт как K →
% Throughput монотонно убывает с ростом K. Это физически корректно.
Thr_e2e_hops = (1 - PER_single).^hop_range .* payloadBits ./ (hop_range * T_pkt_s);

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
        % FastMode=false — честный Doppler-sweep
        modem_d = LoRaModem(fc_Hz, sf, bw, fs, ...
            'CR', CR, 'HasHeader', true, 'UseCRC', true, ...
            'PreambleLen', 8, 'FastMode', false);
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
% Показываем оба канала: AWGN (теоретический минимум) и Rayleigh TDL (городской).
% Важно: для Rayleigh-канала PER насыщается на уровне ~0.1-0.3 при высоком SNR —
% это физически корректно и называется "diversity floor" (предел разнесения).
% Причина: один случайный коэффициент h на пакет (медленные замирания) →
% часть пакетов попадает в глубокие ямы |h|^2 << 1 независимо от SNR.
figure('Name','PER vs SNR','Color','w','Position',[50 50 680 480]);
semilogy(snr_sweep, max(PER_awgn, 1e-4), '-o', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 6);
hold on;
semilogy(snr_sweep, max(PER_tdl,  1e-4), '-^', 'Color', colors(3,:), 'LineWidth', 2.5, 'MarkerSize', 7);
xline(snrThreshold, '--k', 'LineWidth', 1.3, ...
    'Label', sprintf('SF%d threshold = %.1f dB', sf, snrThreshold), ...
    'LabelHorizontalAlignment', 'left');
xlabel('SNR, dB', 'FontSize', 13); ylabel('PER', 'FontSize', 13);
title(sprintf('PER vs SNR | SF=%d, BW=%d kHz, CR=4/%d', sf, bw/1e3, CR+4), 'FontSize', 13);
legend('AWGN (ideal)', 'Rayleigh TDL (urban, diversity floor)', 'Location', 'southwest');
grid on; ylim([1e-4 1]); xlim([snr_sweep(1) snr_sweep(end)]); hold off;

% --- График 2: BER vs SNR ---
figure('Name','BER vs SNR','Color','w','Position',[70 70 680 480]);
semilogy(snr_sweep, max(BER_awgn, 1e-5), '-o', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 6);
hold on;
semilogy(snr_sweep, max(BER_tdl,  1e-5), '-^', 'Color', colors(3,:), 'LineWidth', 2.5, 'MarkerSize', 7);
xline(snrThreshold, '--k', 'LineWidth', 1.3, ...
    'Label', sprintf('SF%d threshold = %.1f dB', sf, snrThreshold), ...
    'LabelHorizontalAlignment', 'left');
xlabel('SNR, dB', 'FontSize', 13); ylabel('BER', 'FontSize', 13);
title(sprintf('BER vs SNR | SF=%d, BW=%d kHz, CR=4/%d', sf, bw/1e3, CR+4), 'FontSize', 13);
legend('AWGN (ideal)', 'Rayleigh TDL (urban, diversity floor)', 'Location', 'southwest');
grid on; ylim([1e-5 1]); xlim([snr_sweep(1) snr_sweep(end)]); hold off;

% --- График 3а: SNR между узлами маршрута (стиль bar, подписи N→M) ---
% Условие nHops > 1 убрано — график строится всегда, даже при одном хопе.
% YDir reverse убран: при отрицательных SNR инверсия уводила столбцы
% за пределы видимой области.
tickLbls = arrayfun(@(h) sprintf('%d\x2192%d', route(h), route(h+1)), ...
    1:nHops, 'UniformOutput', false);

figure('Name', 'SNR between route nodes', 'Color', 'w', 'Position', [90 90 680 460]);
b = bar(1:nHops, routeSNR, 0.6);
b.FaceColor = colors(1,:);   % синий — совпадает с цветовой схемой проекта
b.EdgeColor = 'none';
hold on;

% Линия порога чувствительности
yline(snrThreshold, '--r', 'LineWidth', 1.4, ...
    'Label', sprintf('SNR threshold (SF%d) = %.1f dB', sf, snrThreshold), ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');

% Подписи значений над/под столбцами в зависимости от знака SNR:
% при положительном SNR — текст над столбцом, при отрицательном — под верхушкой
for h = 1:nHops
    if routeSNR(h) >= 0
        yPos   = routeSNR(h) + 0.4;
        vAlign = 'bottom';
    else
        yPos   = routeSNR(h) - 0.4;
        vAlign = 'top';
    end
    text(h, yPos, sprintf('%.1f dB', routeSNR(h)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   vAlign, ...
        'FontSize', 9, 'Color', [0.15 0.15 0.15]);
end

set(gca, 'XTick', 1:nHops, 'XTickLabel', tickLbls, ...
    'XTickLabelRotation', 0, 'FontSize', 11, 'Box', 'on');
xlabel('Hops', 'FontSize', 13);
ylabel('SNR, dB', 'FontSize', 13);
title('SNR between route nodes', 'FontSize', 13, 'FontWeight', 'bold');

% Диапазон оси Y: запас сверху и снизу, порог всегда виден
y_lo = min(min(routeSNR) - 3, snrThreshold - 2);
y_hi = max(routeSNR) * 1.15 + 1;
if y_hi <= y_lo, y_hi = y_lo + 10; end
ylim([y_lo, y_hi]);

grid on; grid minor; hold off;

% --- График 3б: PER vs SNR по хопам маршрута ---
figure('Name', 'PER vs Hop SNR', 'Color', 'w', 'Position', [110 110 660 460]);
scatter(routeSNR, hopPER, 100, colors(2,:), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
hold on;
for h = 1:nHops
    text(routeSNR(h) + 0.25, hopPER(h) + 0.015, ...
        sprintf('%d\x2192%d', route(h), route(h+1)), ...
        'FontSize', 10, 'Color', [0.2 0.2 0.2]);
end
xline(snrThreshold, '--r', 'LineWidth', 1.3, ...
    'Label', sprintf('SNR threshold (SF%d)', sf), ...
    'LabelHorizontalAlignment', 'left');
xlabel('Hop SNR, dB', 'FontSize', 13);
ylabel('Hop PER', 'FontSize', 13);
title(sprintf('PER vs Hop SNR | Route %d\x2192%d | PER_{e2e} = %.4f', ...
    src, dst, PER_e2e), 'FontSize', 13);
grid on;
yMax = max(hopPER) * 1.4 + 0.05;
if isnan(yMax) || yMax < 0.05, yMax = 0.15; end
ylim([0, yMax]); hold off;

% --- График 4: E2E метрики vs число хопов ---
% Цвета осей (YAxis.Color) явно совпадают с цветами соответствующих линий:
%   левая ось  (PER)        — оранжевый colors(2,:)
%   правая ось (Throughput) — синий     colors(1,:)
% Это устраняет путаницу, когда MATLAB автоматически красит подпись оси
% в цвет, не совпадающий с цветом линии на этой оси.
figure('Name', 'E2E Metrics vs Hops', 'Color', 'w', 'Position', [130 130 720 450]);

ax4 = gca;

yyaxis left;
plot(hop_range, PER_e2e_hops, '-o', ...
    'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', colors(2,:));
ylabel('PER_{e2e}', 'FontSize', 13);
ylim([0, 1]);
set(gca, 'YDir', 'normal');
ax4.YAxis(1).Color = colors(2,:);   % левая ось — оранжевая, как линия PER

yyaxis right;
plot(hop_range, Thr_e2e_hops, '-s', ...
    'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', colors(1,:));
ylabel('Throughput_{e2e}, bit/s', 'FontSize', 13);
ylim([0, max(Thr_e2e_hops) * 1.1]);
set(gca, 'YDir', 'normal');
ax4.YAxis(2).Color = colors(1,:);   % правая ось — синяя, как линия Throughput

xline(nHops, '--k', 'LineWidth', 1.4, ...
    'Label', sprintf('Route (%d hops)', nHops), ...
    'LabelVerticalAlignment', 'bottom');
xlabel('Number of hops K', 'FontSize', 13);
title(sprintf('PER_{e2e} and Throughput vs Hops K | worst-case SNR_{hop} = %.1f dB', ...
    snr_typical), 'FontSize', 13);
legend('PER_{e2e} (left axis)', 'Throughput_{e2e} (right axis)', ...
    'Location', 'east', 'FontSize', 11);
grid on;
xlim([1, hop_range(end)]);
set(gca, 'XTick', hop_range);

% --- График 5: PER vs SNR при разных скоростях ---
figure('Name','PER vs SNR — Doppler','Color','w','Position',[130 130 680 480]);
vel_labels = cell(numel(velocities), 1);
for vi = 1:numel(velocities)
    fd_v = (velocities(vi) / 3e8) * fc_Hz;
    semilogy(snr_sweep, max(PER_dynamic(vi,:), 1e-4), '-', 'Color', colors(vi,:), 'LineWidth', 2); hold on;
    vel_labels{vi} = sprintf('v = %d m/s  (f_D = %.0f Hz)', velocities(vi), fd_v);
end
xlabel('SNR, dB', 'FontSize', 13); ylabel('PER', 'FontSize', 13);
title(sprintf('PER vs SNR | SF=%d — Effect of Node Velocity', sf), 'FontSize', 13);
legend(vel_labels, 'Location', 'southwest');
grid on; ylim([1e-3 1]); xlim([snr_sweep(1) snr_sweep(end)]);

% --- График 6: 3D-сеть ---
% Проблема плоского вида: при areaSize=1200 м и высотах Z=5..25 м
% соотношение осей составляет X:Z ≈ 60:1, что делает ось высоты
% практически невидимой при стандартном view(45,30).
% Решение: pbaspect задаёт одинаковый физический масштаб осей X и Y,
% а ось Z масштабируется отдельно (×20) для наглядности высоты.
figure('Name','3D Network','Color','w','Position',[150 150 860 620]);
hold on;

% Серые связи между связными узлами
for i = 1:numNodes
    for j = i+1:numNodes
        if connectivity(i,j)
            plot3([X(i) X(j)],[Y(i) Y(j)],[Z(i) Z(j)], ...
                'Color',[0.80 0.80 0.80],'LineWidth',0.4,'HandleVisibility','off');
        end
    end
end

% Маршрут: цвет хопа кодирует PER (зелёный=0, красный=1)
for h = 1:nHops
    hop_color = [hopPER(h), 1 - hopPER(h), 0];
    plot3(X(route(h:h+1)), Y(route(h:h+1)), Z(route(h:h+1)), ...
        '-', 'Color', hop_color, 'LineWidth', 5);
end

% Узлы: размер и цвет по высоте
scatter3(X, Y, Z, 80, Z, 'filled', 'MarkerEdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.5);

% Выделить src и dst отдельно
plot3(X(src), Y(src), Z(src), 'p', 'MarkerSize', 16, ...
    'MarkerFaceColor', [0 0.6 0], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
plot3(X(dst), Y(dst), Z(dst), 'h', 'MarkerSize', 16, ...
    'MarkerFaceColor', [0.8 0 0], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
text(X(src)+20, Y(src)+20, Z(src)+1, sprintf('src=%d', src), 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0 0.5 0]);
text(X(dst)+20, Y(dst)+20, Z(dst)+1, sprintf('dst=%d', dst), 'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.7 0 0]);

colormap(parula); cb = colorbar;
cb.Label.String = 'Node height, m'; cb.Label.FontSize = 11;

xlabel('X, m', 'FontSize', 12); ylabel('Y, m', 'FontSize', 12); zlabel('Height, m', 'FontSize', 12);
title(sprintf('3D Network | Route %d→%d | %d hops | hop color = PER', src, dst, nHops), 'FontSize', 13);

% Ключевое: pbaspect делает X и Y одинакового масштаба,
% ось Z растягиваем в 20 раз для видимости (иначе 20м против 1200м = невидно)
pbaspect([1, 1, 1/20]);
view(35, 25);    % угол чуть ниже для лучшей видимости Z
grid on; box on; rotate3d on; hold off;

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