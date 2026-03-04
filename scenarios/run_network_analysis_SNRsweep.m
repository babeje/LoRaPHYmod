%% run_network_analysis.m
% Сценарий: оценка энергоэффективности и пропускной способности
% децентрализованной LoRa-сети с учётом динамически меняющихся условий.
%
% Структура сценария:
%   Блок 1  — Параметры сети и топология (базовый SF для связности)
%   Блок 2  — PHY-параметры LoRa + sweep-параметры по SF
%   Блок 3  — SNR-матрица по ITU-R P.1411
%   Блок 4  — Маршрутизация и SNR по хопам
%   Блок 5  — PHY-симуляция: BER/PER по каждому хопу маршрута (SF=7)
%   Блок 6  — End-to-end метрики: PER_e2e, Throughput_e2e
%   Блок 7  — Sweep BER/PER vs SNR по SF={7,9,12}, каналы: AWGN и Rayleigh TDL
%   Блок 8  — PER_e2e и Throughput vs число хопов по SF={7,9,12} (аналитика)
%   Блок 9  — Графики

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
%   areaSize = 4000 м — городской район 4×4 км.
%   numNodes = 50     — плотность ~3.1 узла/км².
%   topoSeed          — автоматически подбирается для максимального числа
%                       хопов, не превышающего hopMax. Связность вычисляется
%                       по порогу SF=7 (snrThreshold = −7.5 дБ) как наиболее
%                       строгому критерию. Маршрут, найденный при SF=7,
%                       используется для всех SF в Блоках 5–8.

numNodes = 50;
areaSize = 4000;   % м

fc_MHz_pre      = 868;
noiseFigure_pre = 6;
txPower_pre     = 14;
snrThresh_pre   = -7.5;      % порог SF=7 — строжайший, базовая топология
k_B_pre         = 1.38e-23;
noisePower_pre  = 10*log10(k_B_pre * 290 * 125e3 * 1000) + noiseFigure_pre;

hopMax = 7;

fprintf('Подбор seed топологии (целевое число хопов: не более %d)...\n', hopMax);

bestSeed  = -1;
bestNhops = 0;

for seedTry = 1:500
    rng(seedTry);
    posT = rand(numNodes, 2) * areaSize;
    rng(seedTry);
    hgtT = 5 + rand(numNodes, 1) * 20;

    connT = false(numNodes);
    for ii = 1:numNodes
        for jj = 1:numNodes
            if ii ~= jj
                dx = posT(ii,1)-posT(jj,1);
                dy = posT(ii,2)-posT(jj,2);
                dz = hgtT(ii)-hgtT(jj);
                d  = sqrt(dx^2+dy^2+dz^2);
                el = atan2d(dz, sqrt(dx^2+dy^2));
                hTx_e = hgtT(ii) + 1.5;
                hRx_e = hgtT(jj) + 1.5;
                L  = ituP1411_corrected(fc_MHz_pre, d, hTx_e, hRx_e, el);
                connT(ii,jj) = (txPower_pre - L - noisePower_pre) > snrThresh_pre;
            end
        end
    end

    maxHops = 0;
    for s = 1:numNodes
        for d = 1:numNodes
            if s ~= d
                r = findRoute(connT, s, d);
                if numel(r)-1 > maxHops
                    maxHops = numel(r)-1;
                end
            end
        end
    end

    if maxHops <= hopMax && maxHops > bestNhops
        bestNhops = maxHops;
        bestSeed  = seedTry;
        fprintf('  seed=%3d → %d хопов  [обновляем]\n', seedTry, maxHops);
        if bestNhops == hopMax
            fprintf('  Достигнут целевой максимум, поиск завершён.\n');
            break;
        end
    end
end

if bestSeed < 0
    warning('Подходящий seed не найден, используется seed=1.');
    bestSeed  = 1;
    bestNhops = 0;
end

topoSeed = bestSeed;
fprintf('Выбран seed=%d  (длиннейший маршрут: %d хопов)\n\n', topoSeed, bestNhops);

rng(topoSeed);
nodePositions = rand(numNodes, 2) * areaSize;
rng(topoSeed);
nodeHeights = 5 + rand(numNodes, 1) * 20;

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
%  БЛОК 2: PHY-ПАРАМЕТРЫ LoRa + SWEEP ПО SF
%% ============================================================
% Sweep по SF реализован через вектор sf_values. Для каждого SF вычисляются
% индивидуальные пороги чувствительности (snrThreshold), длительность
% символа (T_sym_ms) и пакета (T_pkt_s) — все они напрямую влияют на
% форму кривых PER(SNR) и метрики E2E.
%
% Выбор SF={7,9,12} даёт максимальный диапазон компромисса
% «помехоустойчивость vs скорость передачи»:
%   SF7  — максимальная скорость (~5.5 кбит/с), порог −7.5 дБ;
%   SF9  — промежуточный баланс (~2.1 кбит/с), порог −12.5 дБ;
%   SF12 — максимальная дальность (~0.3 кбит/с), порог −20.0 дБ.

sf_values     = [7, 9, 12];          % исследуемые SF
nSF           = numel(sf_values);

fc_Hz        = 868e6;
fc_MHz       = 868;
txPower_dBm  = 14;
noiseFigure  = 6;
hTx_m        = 1.5;
hRx_m        = 1.5;

bw           = 125e3;
fs           = 1e6;
CR           = 1;
payloadBits  = 128;
Npkts        = 1000;
Npkts_sweep  = 1000;

k_B            = 1.38e-23;
T_K            = 290;
noisePower_dBm = 10*log10(k_B * T_K * bw * 1000) + noiseFigure;

% Пороги чувствительности (дБ) для SF 7..12 по спецификации SX1276/SX1262
sf_thresholds  = [-7.5, -10.0, -12.5, -15.0, -17.5, -20.0];

% Вычисление производных параметров для каждого SF
snrThreshold_vec = zeros(1, nSF);
T_sym_ms_vec     = zeros(1, nSF);
T_pkt_ms_vec     = zeros(1, nSF);
T_pkt_s_vec      = zeros(1, nSF);

for si = 1:nSF
    sf_i = sf_values(si);
    snrThreshold_vec(si) = sf_thresholds(sf_i - 6);
    T_sym_ms_vec(si)     = (2^sf_i / bw) * 1e3;
    Nsym_pkt_i           = 8 + ceil((payloadBits/8) * 2 / (sf_i-2)) * (CR+4);
    T_pkt_ms_vec(si)     = (8 + 4.25 + Nsym_pkt_i) * T_sym_ms_vec(si);
    T_pkt_s_vec(si)      = T_pkt_ms_vec(si) / 1e3;
end

% Базовый SF для маршрутизации и графиков хопов (наиболее строгий порог)
sf_base          = 7;
sf_base_idx      = find(sf_values == sf_base, 1);
snrThreshold     = snrThreshold_vec(sf_base_idx);
T_sym_ms         = T_sym_ms_vec(sf_base_idx);
T_pkt_ms         = T_pkt_ms_vec(sf_base_idx);
T_pkt_s          = T_pkt_s_vec(sf_base_idx);

% TDL-профиль: 3 луча, tau_rms ≈ 0.6 мкс (ITU-R P.1411-10, ~300 м)
tdlDelays = [0,  0.5e-6,  1.5e-6];
tdlGains  = [0,  -6,      -12];

% SNR sweep: диапазон перекрывает зону перехода PER 1→0 для всех SF
snr_sweep = -22:1:5;

fprintf('Параметры PHY:\n');
fprintf('  BW=%.0f кГц, CR=4/%d, payload=%d бит\n', bw/1e3, CR+4, payloadBits);
for si = 1:nSF
    fprintf('  SF%-2d: threshold=%.1f дБ, T_sym=%.2f мс, T_pkt≈%.1f мс\n', ...
        sf_values(si), snrThreshold_vec(si), T_sym_ms_vec(si), T_pkt_ms_vec(si));
end
fprintf('  P_noise=%.1f дБм\n\n', noisePower_dBm);

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

% Матрица связности по порогу базового SF=7
connectivity   = snrMatrix > snrThreshold;
connectedNodes = sum(any(connectivity, 2));
snrConn        = snrMatrix(connectivity);

fprintf('Топология (SF=%d, порог=%.1f дБ):\n', sf_base, snrThreshold);
fprintf('  Связанных узлов: %d из %d\n', connectedNodes, numNodes);
fprintf('  SNR: min=%.1f  median=%.1f  max=%.1f дБ\n\n', ...
    min(snrConn), median(snrConn), max(snrConn));

%% ============================================================
%  БЛОК 4: МАРШРУТИЗАЦИЯ И SNR ПО ХОПАМ
%% ============================================================

bestLen = 0;
src = 1; dst = 2; route = [];
for s = 1:numNodes
    for d = 1:numNodes
        if s ~= d
            r = findRoute(connectivity, s, d);
            if numel(r)-1 > bestLen
                bestLen = numel(r)-1;
                src = s; dst = d; route = r;
            end
        end
    end
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
%  БЛОК 5: PHY-СИМУЛЯЦИЯ ПО ХОПАМ МАРШРУТА (базовый SF=7)
%% ============================================================

fprintf('PHY-симуляция по хопам (SF=%d, Npkts=%d)...\n', sf_base, Npkts);

hopBER     = nan(nHops, 1);
hopPER     = nan(nHops, 1);
hopThr_bps = nan(nHops, 1);

for h = 1:nHops
    modem = LoRaModem(fc_Hz, sf_base, bw, fs, ...
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
%  БЛОК 6: END-TO-END МЕТРИКИ (базовый SF=7)
%% ============================================================

PER_e2e         = 1 - prod(1 - hopPER);
Thr_e2e         = min(hopThr_bps);
delay_no_arq_ms = nHops * T_pkt_ms;
delay_arq_ms    = sum(T_pkt_ms ./ max(1 - hopPER, 1e-9));

fprintf('\nEnd-to-end (SF=%d):\n', sf_base);
fprintf('  PER_e2e          = %.4f  (%.1f%%)\n', PER_e2e, PER_e2e*100);
fprintf('  Throughput_e2e   = %.0f бит/с\n',    Thr_e2e);
fprintf('  Задержка (нет ARQ) = %.1f мс\n',     delay_no_arq_ms);
fprintf('  Задержка (ARQ)     = %.1f мс\n\n',   delay_arq_ms);

%% ============================================================
%  БЛОК 7: BER/PER vs SNR — SWEEP ПО SF={7,9,12}, ДВА КАНАЛА
%% ============================================================
% Для каждого SF строятся кривые в двух каналах:
%   AWGN          — теоретическая нижняя граница;
%   Rayleigh TDL  — городской канал с частотной избирательностью.
%
% Обоснование Seed=[]: каждая точка SNR получает независимые реализации
% замираний, что при Npkts_sweep=1000 обеспечивает сходимость к истинному
% математическому ожиданию PER/BER. Фиксированный seed давал зависимость
% результатов от конкретной выборки.

fprintf('Sweep BER/PER vs SNR по SF=%s (%d точек, Npkts=%d)...\n', ...
    mat2str(sf_values), numel(snr_sweep), Npkts_sweep);

% Хранение результатов: строка — SF, столбец — точка SNR
BER_awgn_sf = nan(nSF, numel(snr_sweep));
PER_awgn_sf = nan(nSF, numel(snr_sweep));
BER_tdl_sf  = nan(nSF, numel(snr_sweep));
PER_tdl_sf  = nan(nSF, numel(snr_sweep));

for si = 1:nSF
    sf_i = sf_values(si);
    fprintf('\n  SF=%d (порог=%.1f дБ):\n', sf_i, snrThreshold_vec(si));

    for snri = 1:numel(snr_sweep)
        snr_i = snr_sweep(snri);

        modem_sw = LoRaModem(fc_Hz, sf_i, bw, fs, ...
            'CR', CR, 'HasHeader', true, 'UseCRC', true, ...
            'PreambleLen', 8, 'FastMode', true);

        % Канал AWGN
        [BER_awgn_sf(si, snri), PER_awgn_sf(si, snri), ~] = ...
            LoRaSimulator(modem_sw, AwgnChannel(fs, snr_i)).run(Npkts_sweep, payloadBits);

        % Канал Rayleigh TDL
        ch_tdl = RayleighTDLChannel(fs, snr_i, 0, ...
            'PathDelays', tdlDelays, 'PathGains', tdlGains, 'Seed', []);
        [BER_tdl_sf(si, snri), PER_tdl_sf(si, snri), ~] = ...
            LoRaSimulator(modem_sw, ch_tdl).run(Npkts_sweep, payloadBits);

        fprintf('    SNR=%5.1f дБ | PER: AWGN=%.2f  TDL=%.2f\n', ...
            snr_i, PER_awgn_sf(si, snri), PER_tdl_sf(si, snri));
    end
end

%% ============================================================
%  БЛОК 8: PER_E2E И THROUGHPUT VS ЧИСЛО ХОПОВ — SWEEP ПО SF
%% ============================================================
% Для каждого SF вычисляется PER_single на основе SNR worst-case хопа
% маршрута. Аналитические кривые E2E строятся по формулам:
%   PER_e2e(K) = 1 − (1 − PER_hop)^K
%   Thr_e2e(K) = (1 − PER_hop)^K · payloadBits / (K · T_pkt_s)
%
% Это раскрывает компромисс: SF12 даёт меньший PER_hop при данном SNR,
% но существенно снижает Throughput из-за большего T_pkt.

fprintf('\nАналитика E2E vs K по SF=%s...\n', mat2str(sf_values));

snr_typical = min(routeSNR);   % worst-case SNR маршрута, не зависит от SF

hop_range        = 1:12;
PER_e2e_hops_sf  = nan(nSF, numel(hop_range));
Thr_e2e_hops_sf  = nan(nSF, numel(hop_range));
PER_single_sf    = nan(1, nSF);

for si = 1:nSF
    sf_i     = sf_values(si);
    Tpkt_i   = T_pkt_s_vec(si);

    modem_e2e = LoRaModem(fc_Hz, sf_i, bw, fs, ...
        'CR', CR, 'HasHeader', true, 'UseCRC', true, ...
        'PreambleLen', 8, 'FastMode', true);
    ch_e2e = RayleighTDLChannel(fs, snr_typical, 0, ...
        'PathDelays', tdlDelays, 'PathGains', tdlGains, 'Seed', []);
    [~, PER_single_sf(si), ~] = LoRaSimulator(modem_e2e, ch_e2e).run(Npkts, payloadBits);

    PER_e2e_hops_sf(si, :) = 1 - (1 - PER_single_sf(si)).^hop_range;
    Thr_e2e_hops_sf(si, :) = (1 - PER_single_sf(si)).^hop_range ...
                              .* payloadBits ./ (hop_range * Tpkt_i);

    fprintf('  SF%-2d: PER_hop@SNR_min=%.1f дБ = %.4f  (T_pkt=%.1f мс)\n', ...
        sf_i, snr_typical, PER_single_sf(si), Tpkt_i*1e3);
end

%% ============================================================
%  БЛОК 9: ГРАФИКИ
%% ============================================================

% Цветовая схема: одна линия = один SF
% Маркер: 'o' — AWGN, '^' — Rayleigh TDL
sf_colors = [0.00 0.45 0.70;   % SF7  — синий
             0.85 0.33 0.10;   % SF9  — оранжевый
             0.47 0.67 0.19];  % SF12 — зелёный

sf_line_styles = {'-o', '-s', '-^'};   % стиль линии для SF7, SF9, SF12

colors = [0.00 0.45 0.70;
          0.85 0.33 0.10;
          0.47 0.67 0.19;
          0.63 0.08 0.18];

% --- График 1: PER vs SNR — AWGN и Rayleigh TDL, sweep по SF ---
% Каждый SF — своей цвет. Сплошная линия = AWGN, пунктир = TDL.
% Вертикальные линии — пороги чувствительности соответствующего SF.
figure('Name','PER vs SNR — SF sweep','Color','w','Position',[50 50 780 520]);
hold on;
lgd_entries = {};
for si = 1:nSF
    sf_i = sf_values(si);
    clr  = sf_colors(si, :);

    semilogy(snr_sweep, max(PER_awgn_sf(si,:), 1e-6), '-o', ...
        'Color', clr, 'LineWidth', 2, 'MarkerSize', 5, ...
        'MarkerFaceColor', clr);
    semilogy(snr_sweep, max(PER_tdl_sf(si,:),  1e-6), '--^', ...
        'Color', clr, 'LineWidth', 2, 'MarkerSize', 5);

    xline(snrThreshold_vec(si), ':', 'Color', clr, 'LineWidth', 1.2, ...
        'Label', sprintf('SF%d thr', sf_i), ...
        'LabelHorizontalAlignment', 'right', ...
        'LabelVerticalAlignment', 'bottom', 'FontSize', 8);

    lgd_entries{end+1} = sprintf('SF%d  AWGN', sf_i); %#ok<SAGROW>
    lgd_entries{end+1} = sprintf('SF%d  Rayleigh TDL', sf_i); %#ok<SAGROW>
end
xlabel('SNR, dB', 'FontSize', 13);
ylabel('PER', 'FontSize', 13);
title(sprintf('PER vs SNR | BW=%d kHz, CR=4/%d | SF sweep', bw/1e3, CR+4), 'FontSize', 13);
legend(lgd_entries, 'Location', 'southwest', 'NumColumns', 2, 'FontSize', 9);
grid on; ylim([1e-6 1]); xlim([snr_sweep(1) snr_sweep(end)]); hold off;

% --- График 2: BER vs SNR — AWGN и Rayleigh TDL, sweep по SF ---
figure('Name','BER vs SNR — SF sweep','Color','w','Position',[70 70 780 520]);
hold on;
lgd_entries = {};
for si = 1:nSF
    sf_i = sf_values(si);
    clr  = sf_colors(si, :);

    semilogy(snr_sweep, max(BER_awgn_sf(si,:), 1e-6), '-o', ...
        'Color', clr, 'LineWidth', 2, 'MarkerSize', 5, ...
        'MarkerFaceColor', clr);
    semilogy(snr_sweep, max(BER_tdl_sf(si,:),  1e-6), '--^', ...
        'Color', clr, 'LineWidth', 2, 'MarkerSize', 5);

    xline(snrThreshold_vec(si), ':', 'Color', clr, 'LineWidth', 1.2, ...
        'Label', sprintf('SF%d thr', sf_i), ...
        'LabelHorizontalAlignment', 'right', ...
        'LabelVerticalAlignment', 'bottom', 'FontSize', 8);

    lgd_entries{end+1} = sprintf('SF%d  AWGN', sf_i); %#ok<SAGROW>
    lgd_entries{end+1} = sprintf('SF%d  Rayleigh TDL', sf_i); %#ok<SAGROW>
end
xlabel('SNR, dB', 'FontSize', 13);
ylabel('BER', 'FontSize', 13);
title(sprintf('BER vs SNR | BW=%d kHz, CR=4/%d | SF sweep', bw/1e3, CR+4), 'FontSize', 13);
legend(lgd_entries, 'Location', 'southwest', 'NumColumns', 2, 'FontSize', 9);
grid on; ylim([1e-6 1]); xlim([snr_sweep(1) snr_sweep(end)]); hold off;

% --- График 3а: SNR между узлами маршрута ---
tickLbls = arrayfun(@(h) sprintf('%d\x2192%d', route(h), route(h+1)), ...
    1:nHops, 'UniformOutput', false);

figure('Name', 'SNR between route nodes', 'Color', 'w', 'Position', [90 90 680 460]);
b = bar(1:nHops, routeSNR, 0.6);
b.FaceColor = colors(1,:);
b.EdgeColor = 'none';
hold on;
yline(snrThreshold, '--r', 'LineWidth', 1.4, ...
    'Label', sprintf('SNR threshold (SF%d) = %.1f dB', sf_base, snrThreshold), ...
    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
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
y_lo = min(min(routeSNR) - 3, snrThreshold - 2);
y_hi = max(routeSNR) * 1.15 + 1;
if y_hi <= y_lo, y_hi = y_lo + 10; end
ylim([y_lo, y_hi]);
grid on; grid minor; hold off;

% --- График 3б: PER vs SNR по хопам маршрута (SF=7) ---
figure('Name', 'PER vs Hop SNR (SF base)', 'Color', 'w', 'Position', [110 110 660 460]);
scatter(routeSNR, hopPER, 100, colors(2,:), 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
hold on;
for h = 1:nHops
    text(routeSNR(h) + 0.25, hopPER(h) + 0.015, ...
        sprintf('%d\x2192%d', route(h), route(h+1)), ...
        'FontSize', 10, 'Color', [0.2 0.2 0.2]);
end
xline(snrThreshold, '--r', 'LineWidth', 1.3, ...
    'Label', sprintf('SNR threshold (SF%d)', sf_base), ...
    'LabelHorizontalAlignment', 'left');
xlabel('Hop SNR, dB', 'FontSize', 13);
ylabel('Hop PER', 'FontSize', 13);
title(sprintf('PER vs Hop SNR | SF=%d | Route %d\x2192%d | PER_{e2e} = %.4f', ...
    sf_base, src, dst, PER_e2e), 'FontSize', 13);
grid on;
yMax = max(hopPER) * 1.4 + 0.05;
if isnan(yMax) || yMax < 0.05, yMax = 0.15; end
ylim([0, yMax]); hold off;

% --- График 4: PER_e2e и Throughput vs число хопов K — sweep по SF ---
% Левая ось: PER_e2e(K) — одинаковый диапазон [0..1] для всех SF.
% Правая ось: Throughput(K) — масштаб определяется SF7 (максимальный Thr).
% Каждый SF — своей цвет; сплошные линии = PER, пунктирные = Throughput.
% Вертикальная линия — фактическое число хопов маршрута.
figure('Name', 'E2E Metrics vs Hops — SF sweep', 'Color', 'w', 'Position', [130 130 800 500]);
ax4 = gca;
hold on;

lgd_per = {};
lgd_thr = {};

for si = 1:nSF
    sf_i = sf_values(si);
    clr  = sf_colors(si, :);

    yyaxis left;
    plot(hop_range, PER_e2e_hops_sf(si, :), '-', ...
        'Color', clr, 'LineWidth', 2.0, sf_line_styles{si}{:});

    yyaxis right;
    plot(hop_range, Thr_e2e_hops_sf(si, :), '--', ...
        'Color', clr, 'LineWidth', 1.5);

    lgd_per{end+1} = sprintf('PER_{e2e}  SF%d', sf_i);  %#ok<SAGROW>
    lgd_thr{end+1} = sprintf('Thr_{e2e}  SF%d', sf_i);  %#ok<SAGROW>
end

yyaxis left;
xline(nHops, '--k', 'LineWidth', 1.4, ...
    'Label', sprintf('Route (%d hops)', nHops), ...
    'LabelVerticalAlignment', 'bottom');
ylabel('PER_{e2e}', 'FontSize', 13);
ylim([0, 1]);
ax4.YAxis(1).Color = [0.15 0.15 0.15];

yyaxis right;
ylabel('Throughput_{e2e}, bit/s', 'FontSize', 13);
thr_max = max(Thr_e2e_hops_sf(:));
if thr_max > 0
    ylim([0, thr_max * 1.15]);
end
ax4.YAxis(2).Color = [0.15 0.15 0.15];

xlabel('Number of hops K', 'FontSize', 13);
title(sprintf('PER_{e2e} and Throughput vs Hops K | SNR_{worst-hop} = %.1f dB | SF sweep', ...
    snr_typical), 'FontSize', 13);
legend([lgd_per lgd_thr], 'Location', 'east', 'NumColumns', 2, 'FontSize', 9);
grid on;
xlim([1, hop_range(end)]);
set(gca, 'XTick', hop_range);
hold off;

% --- График 5: 3D-сеть ---
figure('Name','3D Network','Color','w','Position',[150 150 860 620]);
hold on;
for i = 1:numNodes
    for j = i+1:numNodes
        if connectivity(i,j)
            plot3([X(i) X(j)],[Y(i) Y(j)],[Z(i) Z(j)], ...
                'Color',[0.80 0.80 0.80],'LineWidth',0.4,'HandleVisibility','off');
        end
    end
end
for h = 1:nHops
    hop_color = [hopPER(h), 1 - hopPER(h), 0];
    plot3(X(route(h:h+1)), Y(route(h:h+1)), Z(route(h:h+1)), ...
        '-', 'Color', hop_color, 'LineWidth', 5);
end
scatter3(X, Y, Z, 80, Z, 'filled', 'MarkerEdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.5);
plot3(X(src), Y(src), Z(src), 'p', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0 0.6 0], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
plot3(X(dst), Y(dst), Z(dst), 'h', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0.8 0 0], 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
text(X(src)+20, Y(src)+20, Z(src)+1, sprintf('src=%d', src), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', [0 0.5 0]);
text(X(dst)+20, Y(dst)+20, Z(dst)+1, sprintf('dst=%d', dst), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.7 0 0]);
colormap(parula); cb = colorbar;
cb.Label.String = 'Node height, m'; cb.Label.FontSize = 11;
xlabel('X, m', 'FontSize', 12); ylabel('Y, m', 'FontSize', 12); zlabel('Height, m', 'FontSize', 12);
title(sprintf('3D Network | Route %d→%d | %d hops | hop color = PER (SF%d)', ...
    src, dst, nHops, sf_base), 'FontSize', 13);
pbaspect([1, 1, 1/20]);
view(35, 25);
grid on; box on; rotate3d on; hold off;

%% ============================================================
%  СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
%% ============================================================

outDir = fullfile(projectRoot, 'results', 'data');
if ~exist(outDir, 'dir'), mkdir(outDir); end

ts    = datestr(now, 'yyyymmdd_HHMMSS');
fname = fullfile(outDir, ['network_analysis_' ts '.mat']);
save(fname, ...
    'sf_values', 'snr_sweep', ...
    'BER_awgn_sf', 'BER_tdl_sf', ...
    'PER_awgn_sf', 'PER_tdl_sf', ...
    'routeSNR', 'routeDist', 'hopBER', 'hopPER', 'hopThr_bps', ...
    'PER_e2e', 'Thr_e2e', 'delay_no_arq_ms', 'delay_arq_ms', ...
    'PER_e2e_hops_sf', 'Thr_e2e_hops_sf', 'PER_single_sf', 'hop_range', ...
    'route', 'src', 'dst', 'nHops', ...
    'sf_values', 'bw', 'payloadBits', 'Npkts', 'T_pkt_ms_vec', ...
    'tdlDelays', 'tdlGains', 'snr_typical', ...
    'topoSeed', 'numNodes', 'areaSize');

fprintf('\nРезультаты сохранены: %s\n', fname);
fprintf('=== Симуляция завершена ===\n');
