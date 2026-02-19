%% run_network_analysis.m
% Сценарий: оценка энергоэффективности и пропускной способности
% децентрализованной LoRa-сети с учётом динамически меняющихся условий.
%
% Структура сценария:
%   Блок 1 — Параметры сети и топология (совместимо с работой Цоя и др.)
%   Блок 2 — PHY-параметры LoRa и настройки симуляции
%   Блок 3 — Построение SNR-матрицы (ITU-R P.1411, корректные параметры)
%   Блок 4 — Маршрутизация и извлечение SNR по хопам
%   Блок 5 — PHY-симуляция: BER/PER для каждого хопа маршрута
%   Блок 6 — End-to-end метрики: PER_e2e, Throughput_e2e
%   Блок 7 — Кривые BER/PER vs SNR (один хоп, три канала)
%   Блок 8 — Влияние числа хопов на PER_e2e и Throughput
%   Блок 9 — Влияние динамики канала (v=0 vs v=30 м/с)
%   Блок 10 — Графики

clear; clc; close all;

% Добавляем пути к классам проекта
projectRoot = fileparts(mfilename('fullpath'));
projectRoot = fileparts(projectRoot);
addpath(genpath(projectRoot));
addpath(genpath(fullfile(projectRoot, '..', 'itupr1411light_fixpoint')));

fprintf('=== Анализ LoRa mesh-сети: PHY-уровень ===\n\n');

%% ============================================================
%  БЛОК 1: ПАРАМЕТРЫ СЕТИ И ТОПОЛОГИЯ
%% ============================================================

numNodes = 50;
areaSize = 800;          % м

rng(42);
nodePositions = rand(numNodes, 2) * areaSize;

% Высоты узлов (метры над землёй — физически корректно для ITU-R P.1411)
nodeHeights = [38 141 73 138 134 80 145 130 72 135 26 21 102 147 140 143 51 69 97 ...
    87 142 138 8 64 62 23 58 2 129 61 9 142 116 122 31 72 132 150 50 58 4 25 44 77 27 53 81 110 116 42]';

X = nodePositions(:, 1);
Y = nodePositions(:, 2);
Z = nodeHeights;

% 3D-расстояния и углы места между всеми парами узлов
distances  = zeros(numNodes);
elevAngles = zeros(numNodes);
for i = 1:numNodes
    for j = 1:numNodes
        dx = X(i) - X(j);
        dy = Y(i) - Y(j);
        dz = Z(i) - Z(j);
        distances(i, j)  = sqrt(dx^2 + dy^2 + dz^2);
        elevAngles(i, j) = atan2d(dz, sqrt(dx^2 + dy^2));
    end
end

%% ============================================================
%  БЛОК 2: PHY-ПАРАМЕТРЫ LoRa
%% ============================================================

% --- Параметры радиоканала ---
fc_MHz      = 868;          % несущая частота, МГц (стандарт LoRa EU868)
fc_Hz       = fc_MHz * 1e6; % несущая частота, Гц
txPower_dBm = 14;           % типовая мощность TX для EU868, дБм
noiseFigure = 6;            % шум-фигура приёмника, дБ (типовой для SX1276)
hTx_m       = 1.5;          % высота антенны TX над уровнем узла, м
hRx_m       = 1.5;          % высота антенны RX над уровнем узла, м
snrThreshold = -7.5;        % порог чувствительности LoRa SF7, дБ

% --- Параметры LoRa PHY ---
sf          = 7;
bw          = 125e3;        % Гц
fs          = 1e6;          % частота дискретизации, Гц
CR          = 1;            % code rate 4/5
payloadBits = 64;           % длина полезной нагрузки, бит
Npkts       = 50;           % пакетов на точку симуляции

% --- Расчёт тепловой мощности шума ---
k_B        = 1.38e-23;     % постоянная Больцмана
T_K        = 290;          % температура, К
noisePower_dBm = 10*log10(k_B * T_K * bw * 1000) + noiseFigure;

% --- Длительность пакета (для расчёта Throughput) ---
% LoRa: T_pkt = (Npreamble + 4.25 + Nsym_payload) * T_sym
% Оценка: для SF7, BW=125кГц, 8 байт payload + header
T_sym_ms    = (2^sf / bw) * 1e3;           % мс, длительность символа
Nsym_pkt    = 8 + ceil((payloadBits/8) * 2 / (sf - 2)) * (CR + 4);
T_pkt_ms    = (8 + 4.25 + Nsym_pkt) * T_sym_ms;  % полное время пакета, мс
T_pkt_s     = T_pkt_ms / 1e3;

fprintf('Параметры PHY:\n');
fprintf('  SF=%d, BW=%.0f кГц, CR=4/%d\n', sf, bw/1e3, CR+4);
fprintf('  Длительность символа: %.3f мс\n', T_sym_ms);
fprintf('  Длительность пакета:  ~%.1f мс\n', T_pkt_ms);
fprintf('  Порог чувствительности: %.1f дБ\n', snrThreshold);
fprintf('  Мощность шума: %.1f дБм\n\n', noisePower_dBm);

%% ============================================================
%  БЛОК 3: SNR-МАТРИЦА (ITU-R P.1411, корректные параметры)
%% ============================================================

pathLoss  = zeros(numNodes);
snrMatrix = zeros(numNodes);

for i = 1:numNodes
    for j = 1:numNodes
        if i ~= j
            % Эффективная высота антенны = высота узла + высота подвеса антенны
            hTx_eff = Z(i) + hTx_m;
            hRx_eff = Z(j) + hRx_m;

            L = ituP1411_corrected(fc_MHz, distances(i, j), ...
                hTx_eff, hRx_eff, elevAngles(i, j));

            pathLoss(i, j)  = L;
            rxPower         = txPower_dBm - L;
            snrMatrix(i, j) = rxPower - noisePower_dBm;
        else
            snrMatrix(i, j) = NaN;
        end
    end
end

connectivity   = snrMatrix > snrThreshold;
connectedNodes = sum(any(connectivity, 2));
fprintf('Топология сети:\n');
fprintf('  Связанных узлов: %d из %d\n', connectedNodes, numNodes);
fprintf('  SNR min (связные): %.1f дБ\n', min(snrMatrix(connectivity)));
fprintf('  SNR max (связные): %.1f дБ\n\n', max(snrMatrix(connectivity)));

%% ============================================================
%  БЛОК 4: МАРШРУТИЗАЦИЯ И SNR ПО ХОПАМ
%% ============================================================

% Поиск наиболее удалённой пары (совместимо с оригинальным main.m)
[src, dst] = findMostDistantNodes(nodePositions);

% BFS-маршрут по матрице связности
route = findRoute(connectivity, src, dst);

if isempty(route)
    warning('Маршрут между узлами %d и %d не найден.', src, dst);
    route = [];
    nHops = 0;
    routeSNR = [];
else
    nHops     = length(route) - 1;
    routeSNR  = zeros(nHops, 1);
    routeDist = zeros(nHops, 1);

    fprintf('Маршрут src=%d → dst=%d:\n', src, dst);
    for h = 1:nHops
        n1 = route(h);
        n2 = route(h + 1);
        routeSNR(h)  = snrMatrix(n1, n2);
        routeDist(h) = distances(n1, n2);
        fprintf('  Хоп %d: узел %2d → %2d  |  dist=%.0f м  |  SNR=%.1f дБ\n', ...
            h, n1, n2, routeDist(h), routeSNR(h));
    end
    fprintf('  Итого хопов: %d\n\n', nHops);
end

%% ============================================================
%  БЛОК 5: PHY-СИМУЛЯЦИЯ ПО ХОПАМ МАРШРУТА
%% ============================================================
% Для каждого хопа запускаем LoRaSimulator с реальным SNR из матрицы.
% Канал: RayleighTDLChannel (частотно-селективный, медленные замирания).

fprintf('PHY-симуляция по хопам маршрута (Npkts=%d)...\n', Npkts);

hopBER     = nan(nHops, 1);
hopPER     = nan(nHops, 1);
hopThr_bps = nan(nHops, 1);   % пропускная способность хопа, бит/с

for h = 1:nHops
    snr_h = routeSNR(h);

    modem = LoRaModem(fc_Hz, sf, bw, fs, ...
        'CR', CR, ...
        'HasHeader', true, ...
        'UseCRC', true, ...
        'PreambleLen', 8, ...
        'FastMode', true);      % FastMode=true для скорости sweep

    channel = RayleighTDLChannel(fs, snr_h, 0, ...
        'PathDelays', [0, 1.0e-6, 4.0e-6], ...
        'PathGains',  [0, -3,     -6], ...
        'Seed',       42 + h);  % разный seed для каждого хопа

    sim = LoRaSimulator(modem, channel);
    [ber_h, per_h, ~] = sim.run(Npkts, payloadBits);

    hopBER(h)     = ber_h;
    hopPER(h)     = per_h;
    % Throughput хопа с учётом PER: успешные биты в единицу времени
    hopThr_bps(h) = (1 - per_h) * payloadBits / T_pkt_s;

    fprintf('  Хоп %d: SNR=%.1f дБ  BER=%.2e  PER=%.3f  Thr=%.0f бит/с\n', ...
        h, snr_h, ber_h, per_h, hopThr_bps(h));
end

%% ============================================================
%  БЛОК 6: END-TO-END МЕТРИКИ
%% ============================================================
% PER_e2e: пакет потерян если хотя бы один хоп его потерял.
%   PER_e2e = 1 - prod(1 - PER_hop)
%
% Throughput_e2e: минимальный по узким горлам маршрута.
%   Для цепочки хопов с ARQ на каждом хопе:
%   Thr_e2e = payloadBits / (nHops * T_pkt / (1 - PER_hop_avg))
%   Используем "узкое горло" — минимальный Throughput среди хопов.

if nHops > 0
    PER_e2e  = 1 - prod(1 - hopPER);
    Thr_e2e  = min(hopThr_bps);    % bottleneck

    % Задержка: нижняя оценка без ретрансмиссий
    delay_ms = nHops * T_pkt_ms;
    % С учётом ARQ-ретрансмиссий (среднее число попыток = 1/(1-PER))
    delay_arq_ms = sum(T_pkt_ms ./ (1 - max(hopPER, 1e-6)));

    fprintf('\nEnd-to-end метрики маршрута:\n');
    fprintf('  PER_e2e       = %.4f (%.1f%%)\n', PER_e2e, PER_e2e*100);
    fprintf('  Throughput    = %.0f бит/с\n', Thr_e2e);
    fprintf('  Задержка (без ARQ) = %.1f мс\n', delay_ms);
    fprintf('  Задержка (с ARQ)   = %.1f мс\n', delay_arq_ms);
end

%% ============================================================
%  БЛОК 7: BER/PER vs SNR — ОДИН ХОП, ТРИ КАНАЛА
%% ============================================================
% Базовые кривые: AWGN, плоский Rayleigh, TDL Rayleigh.
% Это ключевой результат, показывающий влияние типа канала.

fprintf('\nСнятие кривых BER/PER vs SNR (три канала)...\n');

snr_sweep = -15:1:5;       % дБ
Npkts_sweep = 30;           % пакетов на точку (увеличить для финала)

BER_awgn     = nan(size(snr_sweep));
BER_rayleigh = nan(size(snr_sweep));
BER_tdl      = nan(size(snr_sweep));
PER_awgn     = nan(size(snr_sweep));
PER_rayleigh = nan(size(snr_sweep));
PER_tdl      = nan(size(snr_sweep));

for si = 1:numel(snr_sweep)
    snr_i = snr_sweep(si);

    modem = LoRaModem(fc_Hz, sf, bw, fs, ...
        'CR', CR, 'HasHeader', true, 'UseCRC', true, ...
        'PreambleLen', 8, 'FastMode', true);

    % --- AWGN ---
    ch_awgn = AwgnChannel(fs, snr_i);
    [BER_awgn(si), PER_awgn(si), ~] = LoRaSimulator(modem, ch_awgn).run(Npkts_sweep, payloadBits);

    % --- Плоский Rayleigh (один луч, нет задержек) ---
    ch_ray = RayleighTDLChannel(fs, snr_i, 0, ...
        'PathDelays', 0, ...
        'PathGains',  0, ...
        'Seed', 42);
    [BER_rayleigh(si), PER_rayleigh(si), ~] = LoRaSimulator(modem, ch_ray).run(Npkts_sweep, payloadBits);

    % --- TDL Rayleigh (городской канал, 3 луча) ---
    ch_tdl = RayleighTDLChannel(fs, snr_i, 0, ...
        'PathDelays', [0, 1.0e-6, 4.0e-6], ...
        'PathGains',  [0, -3,     -6], ...
        'Seed', 42);
    [BER_tdl(si), PER_tdl(si), ~] = LoRaSimulator(modem, ch_tdl).run(Npkts_sweep, payloadBits);

    fprintf('  SNR=%5.1f дБ  |  PER: AWGN=%.2f  Ray=%.2f  TDL=%.2f\n', ...
        snr_i, PER_awgn(si), PER_rayleigh(si), PER_tdl(si));
end

%% ============================================================
%  БЛОК 8: PER_E2E И THROUGHPUT КАК ФУНКЦИЯ ЧИСЛА ХОПОВ
%% ============================================================
% Фиксируем типовой SNR хопа (медиана по маршруту), варьируем K хопов.

fprintf('\nВлияние числа хопов на PER_e2e и Throughput...\n');

if nHops > 0
    snr_typical = median(routeSNR);
else
    snr_typical = 5;   % дБ, запасной вариант
end

hop_range  = 1:10;
PER_e2e_hops = nan(size(hop_range));
Thr_e2e_hops = nan(size(hop_range));

% Симулируем один хоп с типовым SNR, затем масштабируем аналитически
modem_h = LoRaModem(fc_Hz, sf, bw, fs, ...
    'CR', CR, 'HasHeader', true, 'UseCRC', true, ...
    'PreambleLen', 8, 'FastMode', true);
ch_h = RayleighTDLChannel(fs, snr_typical, 0, ...
    'PathDelays', [0, 1.0e-6, 4.0e-6], ...
    'PathGains',  [0, -3,     -6], ...
    'Seed', 42);
[~, PER_single, ~] = LoRaSimulator(modem_h, ch_h).run(Npkts, payloadBits);

for ki = 1:numel(hop_range)
    K = hop_range(ki);
    PER_e2e_hops(ki)  = 1 - (1 - PER_single)^K;
    % Throughput: пропускная способность узкого горла с ARQ
    Thr_e2e_hops(ki)  = (1 - PER_single) * payloadBits / (K * T_pkt_s);
end

fprintf('  SNR типового хопа: %.1f дБ, PER_hop=%.3f\n', snr_typical, PER_single);

%% ============================================================
%  БЛОК 9: ВЛИЯНИЕ ДИНАМИКИ КАНАЛА (v=0 vs v=30 м/с)
%% ============================================================
% Сравниваем статичный и подвижный сценарии по PER.

fprintf('\nВлияние скорости движения на PER...\n');

velocities  = [0, 5, 15, 30];   % м/с
PER_dynamic = nan(numel(velocities), numel(snr_sweep));

for vi = 1:numel(velocities)
    v = velocities(vi);
    fd = (v / 3e8) * fc_Hz;    % допплеровский сдвиг, Гц

    for si = 1:numel(snr_sweep)
        snr_i = snr_sweep(si);

        modem_d = LoRaModem(fc_Hz, sf, bw, fs, ...
            'CR', CR, 'HasHeader', true, 'UseCRC', true, ...
            'PreambleLen', 8, 'FastMode', true);

        % Допплер моделируется через DopplerChannel поверх TDL:
        % применяем DopplerChannel с CFO=0 и заданной скоростью
        ch_d = DopplerChannel(fs, snr_i, 0, fc_Hz, v, 0, 0);

        [~, PER_dynamic(vi, si), ~] = LoRaSimulator(modem_d, ch_d).run(Npkts_sweep, payloadBits);
    end
    fprintf('  v=%2d м/с (fd=%.1f Гц): PER при SNR=0дБ = %.3f\n', ...
        v, fd, PER_dynamic(vi, find(snr_sweep==0, 1)));
end

%% ============================================================
%  БЛОК 10: ГРАФИКИ
%% ============================================================

colors = [0.0 0.45 0.70;   % синий
          0.85 0.33 0.10;  % оранжевый
          0.47 0.67 0.19;  % зелёный
          0.63 0.08 0.18]; % тёмно-красный

% --- График 1: PER vs SNR, три типа канала ---
figure('Name', 'PER vs SNR — типы канала', 'Color', 'w', 'Position', [100 100 680 480]);
semilogy(snr_sweep, max(PER_awgn,     1e-4), '-o', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 6); hold on;
semilogy(snr_sweep, max(PER_rayleigh, 1e-4), '-s', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 6);
semilogy(snr_sweep, max(PER_tdl,      1e-4), '-^', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 6);
xline(snrThreshold, '--k', 'LineWidth', 1.2, 'Label', 'Порог чувствительности');
xlabel('SNR, дБ', 'FontSize', 13);
ylabel('PER', 'FontSize', 13);
title(sprintf('PER vs SNR, SF=%d, BW=%.0f кГц', sf, bw/1e3), 'FontSize', 13);
legend('AWGN', 'Rayleigh (плоский)', 'Rayleigh TDL (городской)', 'Location', 'southwest');
grid on; ylim([1e-3 1]); xlim([snr_sweep(1) snr_sweep(end)]);

% --- График 2: BER vs SNR, три типа канала ---
figure('Name', 'BER vs SNR — типы канала', 'Color', 'w', 'Position', [120 120 680 480]);
semilogy(snr_sweep, max(BER_awgn,     1e-5), '-o', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 6); hold on;
semilogy(snr_sweep, max(BER_rayleigh, 1e-5), '-s', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 6);
semilogy(snr_sweep, max(BER_tdl,      1e-5), '-^', 'Color', colors(3,:), 'LineWidth', 2, 'MarkerSize', 6);
xlabel('SNR, дБ', 'FontSize', 13);
ylabel('BER', 'FontSize', 13);
title(sprintf('BER vs SNR, SF=%d, BW=%.0f кГц', sf, bw/1e3), 'FontSize', 13);
legend('AWGN', 'Rayleigh (плоский)', 'Rayleigh TDL (городской)', 'Location', 'southwest');
grid on; ylim([1e-4 1]); xlim([snr_sweep(1) snr_sweep(end)]);

% --- График 3: SNR и PER по хопам реального маршрута ---
if nHops > 0
    figure('Name', 'SNR и PER по хопам маршрута', 'Color', 'w', 'Position', [140 140 750 420]);
    yyaxis left;
    bar(1:nHops, routeSNR, 0.5, 'FaceColor', colors(1,:), 'FaceAlpha', 0.7);
    ylabel('SNR хопа, дБ', 'FontSize', 13);
    ylim([min(routeSNR)-5, max(routeSNR)+5]);
    yyaxis right;
    plot(1:nHops, hopPER, '-s', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('PER хопа', 'FontSize', 13);
    ylim([0, 1]);
    xlabel('Номер хопа', 'FontSize', 13);
    title(sprintf('Маршрут: узел %d → %d (%d хопов)', src, dst, nHops), 'FontSize', 13);
    tickLabels = arrayfun(@(h) sprintf('%d→%d', route(h), route(h+1)), 1:nHops, 'UniformOutput', false);
    set(gca, 'XTick', 1:nHops, 'XTickLabel', tickLabels, 'XTickLabelRotation', 30);
    legend('SNR хопа', 'PER хопа', 'Location', 'best');
    grid on;
end

% --- График 4: PER_e2e и Throughput vs число хопов ---
figure('Name', 'E2E метрики vs число хопов', 'Color', 'w', 'Position', [160 160 750 420]);
yyaxis left;
plot(hop_range, PER_e2e_hops, '-o', 'Color', colors(2,:), 'LineWidth', 2, 'MarkerSize', 8);
ylabel('PER_{e2e}', 'FontSize', 13);
ylim([0, 1]);
yyaxis right;
plot(hop_range, Thr_e2e_hops, '-s', 'Color', colors(1,:), 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Пропускная способность, бит/с', 'FontSize', 13);
xlabel('Число хопов K', 'FontSize', 13);
title(sprintf('E2E метрики vs число хопов, SNR_{хоп}=%.1f дБ', snr_typical), 'FontSize', 13);
legend('PER_{e2e}', 'Throughput_{e2e}', 'Location', 'east');
grid on; xlim([1, hop_range(end)]);

% --- График 5: PER vs SNR при разных скоростях ---
figure('Name', 'PER vs SNR — влияние скорости', 'Color', 'w', 'Position', [180 180 680 480]);
vel_labels = cell(numel(velocities), 1);
for vi = 1:numel(velocities)
    fd_v = (velocities(vi) / 3e8) * fc_Hz;
    semilogy(snr_sweep, max(PER_dynamic(vi,:), 1e-4), ...
        '-', 'Color', colors(vi,:), 'LineWidth', 2); hold on;
    vel_labels{vi} = sprintf('v=%d м/с (f_D=%.0f Гц)', velocities(vi), fd_v);
end
xlabel('SNR, дБ', 'FontSize', 13);
ylabel('PER', 'FontSize', 13);
title(sprintf('PER vs SNR, SF=%d — влияние скорости движения', sf), 'FontSize', 13);
legend(vel_labels, 'Location', 'southwest');
grid on; ylim([1e-3 1]); xlim([snr_sweep(1) snr_sweep(end)]);

% --- График 6: 3D-сеть с маршрутом (совместимо с оригиналом) ---
figure('Name', '3D Network Graph', 'Color', 'w', 'Position', [200 200 800 600]);
hold on;
for i = 1:numNodes
    for j = i+1:numNodes
        if connectivity(i, j)
            plot3([X(i) X(j)], [Y(i) Y(j)], [Z(i) Z(j)], ...
                'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
    end
end
if ~isempty(route)
    for h = 1:nHops
        % Окраска хопа по его PER: зелёный=хорошо, красный=плохо
        per_color = [hopPER(h), 1-hopPER(h), 0];
        plot3(X(route(h:h+1)), Y(route(h:h+1)), Z(route(h:h+1)), ...
            '-', 'Color', per_color, 'LineWidth', 4);
    end
end
scatter3(X, Y, Z, 60, Z, 'filled');
colormap(jet); colorbar;
xlabel('X, м'); ylabel('Y, м'); zlabel('Высота, м');
title('3D-сеть: связи, маршрут (цвет = PER хопа)', 'FontSize', 13);
grid on; box on; view(45, 30); rotate3d on;
hold off;

%% ============================================================
%  СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
%% ============================================================

outDir = fullfile(projectRoot, 'results', 'data');
if ~exist(outDir, 'dir'), mkdir(outDir); end

ts = datestr(now, 'yyyymmdd_HHMMSS');
save(fullfile(outDir, ['network_analysis_' ts '.mat']), ...
    'snr_sweep', ...
    'BER_awgn', 'BER_rayleigh', 'BER_tdl', ...
    'PER_awgn', 'PER_rayleigh', 'PER_tdl', ...
    'routeSNR', 'hopBER', 'hopPER', 'hopThr_bps', ...
    'PER_e2e', 'Thr_e2e', ...
    'PER_e2e_hops', 'Thr_e2e_hops', 'hop_range', ...
    'PER_dynamic', 'velocities', ...
    'route', 'src', 'dst', 'nHops', ...
    'sf', 'bw', 'payloadBits', 'Npkts', 'T_pkt_ms');

fprintf('\nРезультаты сохранены: %s\n', outDir);
fprintf('=== Симуляция завершена ===\n');