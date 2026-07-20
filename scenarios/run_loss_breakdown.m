% run_loss_breakdown.m
%
% Назначение:
%   Диагностика механизма потерь пакетов в канале с замираниями Релея.
%   Разделяет потери на два физических механизма:
%
%     (A) Потеря преамбулы/синхронизации (ERASURE) — приёмник не обнаружил
%         преамбулу или не смог декодировать заголовок. Теряется ВЕСЬ пакет.
%         Структурно: demodulate() возвращает rx_ok = false.
%
%     (B) Повреждение данных (CRC FAIL) — преамбула обнаружена, данные
%         извлечены нужной длины, но CRC не сошёлся (ошибки в символах).
%         Структурно: demodulate() возвращает rx_ok = true, crc_ok = false.
%
%   Цель: определить, какой механизм доминирует в каждой зоне ОСШ, и тем
%   самым оценить, способна ли оценка канала (Шаг 1+2) снизить потери.
%   Оценка канала помогает против механизма (B) — стирания символов данных
%   в замираниях, восстанавливаемые декодером со стираниями. Против
%   механизма (A) она бессильна — нечего оценивать, если преамбула утонула.
%
% Методология:
%   Плоские замирания Релея AR(1) генерируются ЯВНО (известен истинный h[n]
%   посимвольно), применяются вручную к пакетному сигналу. Это позволяет
%   коррелировать исход пакета с провалами |h| в преамбуле vs данных.
%   Шум нормируется по мощности НЕзамирающего сигнала (фиксированный шумовой
%   порог), поэтому замирание реально снижает мгновенное ОСШ.
%
% Зависимости:
%   core/LoRaModem.m, core/LoRaSimulator.m, phy/LoRaPHY.m
%
% Выходные данные:
%   results/figures/loss_breakdown_per.fig
%   results/figures/loss_breakdown_fade.fig
%   results/data/loss_breakdown_<timestamp>.mat
%
% Запуск:
%   cd scenarios
%   run_loss_breakdown
% ---------------------------------------------------------------

clearvars; close all;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));

%% ============================================================
%  Блок 1 — Параметры PHY и симуляции
%% ============================================================
rf_freq = 868e6;
sf      = 7;
bw      = 125e3;
fs      = 1e6;

CR          = 4;        % code rate: CR=1 → 4/5 (только детектирование)
HasHeader   = true;
UseCRC      = true;
PreambleLen = 8;

% ВАЖНО: FastMode=false — нужна РЕАЛЬНАЯ детекция преамбулы, чтобы провалы
% |h| на преамбуле приводили к ERASURE. В FastMode детекция преамбулы
% может пропускаться, и механизм (A) не проявится.
FastMode = false;

payloadLenBits = 128;
Npkts          = 300;            % пакетов на точку ОСШ
snr_list       = -12 : 2 : 6;    % диапазон ОСШ, дБ

% Параметры подвижности (для AR(1) замираний)
v_mps  = 30;                     % скорость объекта, м/с
fD_Hz  = (v_mps / 3e8) * rf_freq;
T_sym  = 2^sf / bw;
Ns     = round(2^sf / bw * fs);  % отсчётов на символ
rho    = besselj(0, 2*pi * fD_Hz * T_sym);   % коэффициент AR(1)

fprintf('\n=== Диагностика механизма потерь: преамбула vs данные ===\n');
fprintf('SF=%d, BW=%.0f кГц, CR=4/%d, FastMode=%d\n', sf, bw/1e3, CR+4, FastMode);
fprintf('v=%.0f м/с, fD=%.1f Гц, rho=%.4f, Ns=%d\n', v_mps, fD_Hz, rho, Ns);
fprintf('Преамбула=%d симв, пакетов/точку=%d\n\n', PreambleLen, Npkts);

%% ============================================================
%  Блок 2 — Основной цикл по ОСШ
%% ============================================================
N_snr = numel(snr_list);

% Счётчики исходов: [N_snr x 1]
cnt_erasure = zeros(N_snr, 1);   % механизм (A): потеря преамбулы
cnt_crcfail = zeros(N_snr, 1);   % механизм (B): повреждение данных
cnt_success = zeros(N_snr, 1);

% Статистика провалов |h| по исходам (для подтверждения причинности):
% накапливаем минимальное |h| в преамбуле и в данных для каждого исхода
minH_pre_erasure = [];   % min|h| преамбулы для ERASURE-пакетов
minH_dat_crcfail = [];   % min|h| данных для CRC-FAIL пакетов
minH_pre_success = [];   % min|h| преамбулы для SUCCESS-пакетов
minH_dat_success = [];   % min|h| данных для SUCCESS-пакетов

for si = 1 : N_snr
    snr_dB  = snr_list(si);
    snr_lin = 10^(snr_dB / 10);

    modem = LoRaModem(rf_freq, sf, bw, fs, ...
        'CR', CR, 'HasHeader', HasHeader, 'UseCRC', UseCRC, ...
        'PreambleLen', PreambleLen, 'FastMode', FastMode);

    for pk = 1 : Npkts
        % --- Генерация и модуляция пакета ---
        bits_tx = logical(randi([0 1], payloadLenBits, 1));
        [txSig, ~, ~] = modem.modulate(bits_tx);
        txSig = txSig(:);
        L = length(txSig);

        % --- Генерация плоских замираний AR(1) (известный h[n]) ---
        n_sym = ceil(L / Ns);
        h_sym = ar1_flat_fading(n_sym, rho);   % [n_sym x 1], E{|h|^2}=1

        % Развёртка посимвольного h в поотсчётный вектор
        h_vec = kron(h_sym, ones(Ns, 1));
        h_vec = h_vec(1:L);

        % --- Применение замираний + АБГШ ---
        % Шум нормируется по мощности НЕзамирающего сигнала (фиксированный
        % порог) → замирание реально снижает мгновенное ОСШ.
        P_tx    = mean(abs(txSig).^2);
        P_noise = P_tx / snr_lin;
        noise   = sqrt(P_noise/2) * (randn(L,1) + 1j*randn(L,1));
        rxSig   = h_vec .* txSig + noise;

        % --- Демодуляция через реальный тракт ---
        [~, rx_ok, crc_ok] = modem.demodulate(rxSig);

        % --- Статистика провалов |h| по областям пакета ---
        % Преамбула: первые PreambleLen символов
        idx_pre = 1 : min(PreambleLen, n_sym);
        % Данные: всё после преамбулы (включает SFD/заголовок — приближённо)
        idx_dat = (PreambleLen+1) : n_sym;

        minH_pre = min(abs(h_sym(idx_pre)));
        if ~isempty(idx_dat)
            minH_dat = min(abs(h_sym(idx_dat)));
        else
            minH_dat = minH_pre;
        end

        % --- Классификация исхода ---
        if ~rx_ok
            % Механизм (A): преамбула/синхронизация/заголовок
            cnt_erasure(si)   = cnt_erasure(si) + 1;
            minH_pre_erasure  = [minH_pre_erasure; minH_pre]; %#ok<AGROW>
        elseif ~crc_ok
            % Механизм (B): данные повреждены
            cnt_crcfail(si)   = cnt_crcfail(si) + 1;
            minH_dat_crcfail  = [minH_dat_crcfail; minH_dat]; %#ok<AGROW>
        else
            cnt_success(si)   = cnt_success(si) + 1;
            minH_pre_success  = [minH_pre_success; minH_pre]; %#ok<AGROW>
            minH_dat_success  = [minH_dat_success; minH_dat]; %#ok<AGROW>
        end
    end

    % --- Сводка по точке ОСШ ---
    per      = (cnt_erasure(si) + cnt_crcfail(si)) / Npkts;
    tot_loss = cnt_erasure(si) + cnt_crcfail(si);
    frac_pre = cnt_erasure(si) / max(tot_loss, 1) * 100;
    frac_dat = cnt_crcfail(si) / max(tot_loss, 1) * 100;

    fprintf('SNR=%4.0f дБ | PER=%.3f | ERASURE=%.3f CRC_FAIL=%.3f | потери: %3.0f%% преамб / %3.0f%% данные\n', ...
        snr_dB, per, cnt_erasure(si)/Npkts, cnt_crcfail(si)/Npkts, frac_pre, frac_dat);
end

fprintf('\n=== Диагностика завершена ===\n\n');

%% ============================================================
%  Блок 3 — Производные метрики
%% ============================================================
per_total   = (cnt_erasure + cnt_crcfail) / Npkts;
per_erasure = cnt_erasure / Npkts;
per_crcfail = cnt_crcfail / Npkts;

%% ============================================================
%  Блок 4 — Сохранение результатов
%% ============================================================
outDir = fullfile(projectRoot, 'results', 'data');
if ~exist(outDir, 'dir'), mkdir(outDir); end
ts = datestr(now, 'yyyymmdd_HHMMSS');
save(fullfile(outDir, sprintf('loss_breakdown_%s.mat', ts)), ...
    'snr_list', 'cnt_erasure', 'cnt_crcfail', 'cnt_success', ...
    'per_total', 'per_erasure', 'per_crcfail', 'Npkts', ...
    'sf', 'bw', 'fs', 'CR', 'v_mps', 'fD_Hz', 'rho', ...
    'minH_pre_erasure', 'minH_dat_crcfail', ...
    'minH_pre_success', 'minH_dat_success');
fprintf('Результаты сохранены: results/data/loss_breakdown_%s.mat\n', ts);

%% ============================================================
%  Блок 5 — График 1: разбивка PER по механизмам
%% ============================================================
figDir = fullfile(projectRoot, 'results', 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

fig1 = figure('Name', 'Разбивка потерь по механизмам', 'Color', 'w', ...
    'Position', [60 60 760 520], 'NumberTitle', 'off');
hold on; grid on; box on;

% Столбчатая диаграмма с накоплением: ERASURE снизу, CRC_FAIL сверху
bar_data = [per_erasure, per_crcfail];
hb = bar(snr_list, bar_data, 0.7, 'stacked');
hb(1).FaceColor = [0.85 0.33 0.10];   % ERASURE — оранжевый
hb(2).FaceColor = [0.00 0.45 0.70];   % CRC_FAIL — синий

% Линия суммарного PER
plot(snr_list, per_total, 'k--o', 'LineWidth', 1.8, 'MarkerSize', 6, ...
    'MarkerFaceColor', 'k');

xlabel('ОСШ (SNR), дБ', 'FontSize', 12);
ylabel('Вероятность исхода', 'FontSize', 12);
title({'Разбивка потерь пакетов по механизмам', ...
    sprintf('SF=%d, BW=%.0f кГц, CR=4/%d, v=%.0f м/с', sf, bw/1e3, CR+4, v_mps)}, ...
    'FontSize', 12);
legend({'ERASURE (потеря преамбулы)', 'CRC FAIL (повреждение данных)', ...
    'PER суммарный'}, 'Location', 'northeast', 'FontSize', 10);
ylim([0 1]);
hold off;
savefig(fig1, fullfile(figDir, 'loss_breakdown_per.fig'));
fprintf('График 1 сохранён: loss_breakdown_per.fig\n');

%% ============================================================
%  Блок 6 — График 2: связь исхода с провалами |h|
%% ============================================================
% Подтверждение причинности: распределение min|h| для каждого исхода.
% Ожидаем: ERASURE-пакеты имеют глубокие провалы в ПРЕАМБУЛЕ,
%          CRC_FAIL-пакеты — глубокие провалы в ДАННЫХ.
fig2 = figure('Name', 'Связь потерь с глубиной замираний', 'Color', 'w', ...
    'Position', [80 80 760 520], 'NumberTitle', 'off');
hold on; grid on; box on;

edges = 0 : 0.1 : 2.0;
if ~isempty(minH_pre_erasure)
    histogram(minH_pre_erasure, edges, 'Normalization', 'probability', ...
        'FaceColor', [0.85 0.33 0.10], 'FaceAlpha', 0.6, ...
        'DisplayName', 'min|h| преамбулы (ERASURE)');
end
if ~isempty(minH_dat_crcfail)
    histogram(minH_dat_crcfail, edges, 'Normalization', 'probability', ...
        'FaceColor', [0.00 0.45 0.70], 'FaceAlpha', 0.6, ...
        'DisplayName', 'min|h| данных (CRC FAIL)');
end
if ~isempty(minH_dat_success)
    histogram(minH_dat_success, edges, 'Normalization', 'probability', ...
        'FaceColor', [0.47 0.67 0.19], 'FaceAlpha', 0.5, ...
        'DisplayName', 'min|h| данных (SUCCESS)');
end

xlabel('Минимальное |h| в области пакета', 'FontSize', 12);
ylabel('Доля пакетов', 'FontSize', 12);
title({'Связь исхода пакета с глубиной замираний', ...
    'Глубокие провалы |h| → потеря в соответствующей области'}, 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 10);
hold off;
savefig(fig2, fullfile(figDir, 'loss_breakdown_fade.fig'));
fprintf('График 2 сохранён: loss_breakdown_fade.fig\n');

%% ============================================================
%  Блок 7 — Итоговая сводка
%% ============================================================
fprintf('\n--- Итоговая разбивка потерь ---\n');
fprintf('%-8s %-8s %-12s %-12s %-12s\n', 'SNR,дБ', 'PER', '%преамбула', '%данные', 'min|h| данных(CRC)');
fprintf('%s\n', repmat('-', 1, 58));
for si = 1 : N_snr
    tot_loss = cnt_erasure(si) + cnt_crcfail(si);
    fp = cnt_erasure(si) / max(tot_loss,1) * 100;
    fd = cnt_crcfail(si) / max(tot_loss,1) * 100;
    fprintf('%-8.0f %-8.3f %-12.0f %-12.0f\n', snr_list(si), per_total(si), fp, fd);
end

if ~isempty(minH_dat_crcfail)
    fprintf('\nМедианный min|h| в данных для CRC-FAIL пакетов: %.3f\n', median(minH_dat_crcfail));
end
if ~isempty(minH_dat_success)
    fprintf('Медианный min|h| в данных для SUCCESS пакетов:  %.3f\n', median(minH_dat_success));
end
fprintf('\nИнтерпретация:\n');
fprintf('  Если в рабочей зоне ОСШ доминирует %% данные → оценка канала\n');
fprintf('  (Шаг 1+2) + erasure-декодирование при CR>=3 могут снизить floor.\n');
fprintf('  Если доминирует %% преамбула → нужно разнесение/повтор.\n\n');

%% ============================================================
%  ЛОКАЛЬНАЯ ФУНКЦИЯ — генератор плоских замираний AR(1)
%% ============================================================
function h = ar1_flat_fading(N_sym, rho)
    % ar1_flat_fading  Плоские замирания Релея по модели AR(1).
    %
    %   h[n] = rho·h[n-1] + sqrt(1-rho^2)·v[n],  v[n] ~ CN(0,1)
    %   h[0] ~ CN(0,1)  (стационарное начальное состояние)
    %   E{|h[n]|^2} = 1  (нормировка мощности)
    %
    % Соответствует модели RayleighTDLChannel в режиме PerSymbol для
    % одного луча (плоский канал, Bc >> BW).

    h = zeros(N_sym, 1);
    h(1) = (randn + 1j*randn) / sqrt(2);
    coeff = sqrt(1 - rho^2);
    for n = 2 : N_sym
        w    = (randn + 1j*randn) / sqrt(2);
        h(n) = rho * h(n-1) + coeff * w;
    end
end
