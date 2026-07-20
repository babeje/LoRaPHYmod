function fig_PER_channel_comparison(varargin)
% fig_PER_channel_comparison  PER vs ОСШ для трёх моделей канала
%
% Каналы: АБГШ / замирания Релея (плоский) / замирания Релея (TDL, 4 луча)
% Параметры по умолчанию: SF=7, BW=125 кГц, CR=4/5, PL=10 байт
%
% Использование:
%   fig_PER_channel_comparison()              % параметры по умолчанию
%   fig_PER_channel_comparison('N_pkts', 500) % задать число пакетов
%
% Аргументы (пары 'имя', значение):
%   'SF'      — коэффициент расширения спектра    (по умолч. 7)
%   'N_pkts'  — число пакетов на точку ОСШ       (по умолч. 1000)
%   'PL'      — длина полезной нагрузки, байт     (по умолч. 10)
%   'SNR_dB'  — вектор значений ОСШ              (по умолч. -14:1:5)
%
% Выход:
%   fig_PER_channel_comparison.png  (300 DPI, академический стиль)

%% ── Разбор аргументов ────────────────────────────────────────────────────
p = inputParser();
addParameter(p, 'SF',     7,          @isnumeric);
addParameter(p, 'N_pkts', 1000,       @isnumeric);
addParameter(p, 'PL',     10,         @isnumeric);
addParameter(p, 'SNR_dB', -14:1:5,   @isnumeric);
parse(p, varargin{:});
SF     = p.Results.SF;
N_pkts = p.Results.N_pkts;
PL     = p.Results.PL;
SNR_dB = p.Results.SNR_dB(:);

%% ── Параметры LoRa ───────────────────────────────────────────────────────
BW  = 125e3;          % полоса модуляции, Гц
M   = 2 ^ SF;         % число символьных состояний
n   = (0 : M-1)';     % индексы отсчётов

% Число символов данных в пакете (CR=4/5)
n_coded  = ceil(PL * 8 * 5/4 / SF) * SF;
N_data   = n_coded / SF;

% Пороговое отношение для сравнения (рисунок)
SNR_thr_dB = -7.5;

%% ── Параметры TDL-канала (ITU-R P.1411-10 Urban, 4 луча) ─────────────────
tau_s   = [0, 0.5e-6, 2.0e-6, 5.0e-6];   % задержки, с
P_dB    = [0, -2,     -5,      -8    ];    % мощности, дБ
P_norm  = 10.^(P_dB/10);
P_norm  = P_norm / sum(P_norm);            % нормировка
sigma_h = sqrt(P_norm / 2);               % СКО вещ./мним. части каждого луча

%% ── CSS модуляция / демодуляция ──────────────────────────────────────────
% Базовый восходящий чирп (s=0):  x_base[n] = exp(j·π·n²/M)
x_base = exp(1j * pi / M * n .^ 2);

% Нисходящий чирп для дечирпирования:
down_chirp = conj(x_base);

% Модуляция символа s: x_s[n] = x_base[n] · exp(j·2π·s·n/M)
modulate = @(s) x_base .* exp(1j * 2 * pi * s / M * n);

% Демодуляция: дечирп → БПФ → argmax (0-индексирован)
demodulate = @(rx) mod(round(M * angle(sum(fft(rx .* down_chirp))) / (2*pi)), M);
% Примечание: используется аргумент суммы, что эквивалентно детектированию
% по индексу максимума |FFT|, оптимизировано для скорости через sum

detect = @(rx) find_peak(fft(rx .* down_chirp));

%% ── Инициализация результатов ────────────────────────────────────────────
PER_awgn = zeros(size(SNR_dB));
PER_rayl = zeros(size(SNR_dB));
PER_tdl  = zeros(size(SNR_dB));

fprintf('SF=%d | BW=%.0f кГц | PL=%d байт | N_data=%d символов | N_pkts=%d\n', ...
        SF, BW/1e3, PL, N_data, N_pkts);
fprintf('%-10s | %-10s | %-15s | %-10s\n', 'SNR, дБ','АБГШ','Релей (плоск.)','Релей TDL');

%% ── Основной цикл симуляции ─────────────────────────────────────────────
for i_snr = 1:length(SNR_dB)

    snr_lin = 10 ^ (SNR_dB(i_snr) / 10);
    % ОСШ в полосе BW: SNR_BW = P_s/(N0*BW)
    % Мощность символа = 1 → СКО шума на отсчёт:
    sigma_n = sqrt(1 / (2 * snr_lin));

    err_awgn = 0;
    err_rayl = 0;
    err_tdl  = 0;

    for i_pkt = 1:N_pkts

        % Случайный пакет
        sym_tx = randi([0, M-1], N_data, 1);
        lost_awgn = false;
        lost_rayl = false;
        lost_tdl  = false;

        for i_sym = 1:N_data

            tx    = modulate(sym_tx(i_sym));
            noise = sigma_n * (randn(M, 1) + 1j * randn(M, 1));

            % ── АБГШ ──────────────────────────────────────────────────────
            rx = tx + noise;
            if detect(rx) ~= sym_tx(i_sym)
                lost_awgn = true;
            end

            % ── Релей плоский (один коэф. на символ) ─────────────────────
            h_flat = sigma_h(1) * (randn + 1j*randn) * sqrt(2);  % CN(0,1)
            rx = h_flat * tx + noise;
            if detect(rx) ~= sym_tx(i_sym)
                lost_rayl = true;
            end

            % ── Релей TDL (4 луча, квазиплоский режим SF=7) ───────────────
            % При τ_max/T_sym ≈ 1.7e-3 (<<1) лучи неразрешимы;
            % результирующий коэффициент = сумма независимых CN-переменных
            h_taps = arrayfun(@(s) s*(randn+1j*randn), sigma_h).';
            % Фазовые сдвиги из-за задержек (на несущей частоте символа):
            % φ_l = 2π·f_c·τ_l, моделируем случайную фазу между символами
            % (для SF=7 f_D·T_sym ≈ 0.09 бина, т.е. изменение между символами мало)
            h_eff = sum(h_taps);  % суперпозиция лучей
            rx = h_eff * tx + noise;
            if detect(rx) ~= sym_tx(i_sym)
                lost_tdl = true;
            end

        end % символы

        err_awgn = err_awgn + lost_awgn;
        err_rayl = err_rayl + lost_rayl;
        err_tdl  = err_tdl  + lost_tdl;

    end % пакеты

    PER_awgn(i_snr) = err_awgn / N_pkts;
    PER_rayl(i_snr) = err_rayl / N_pkts;
    PER_tdl(i_snr)  = err_tdl  / N_pkts;

    fprintf('%-10.1f | %-10.4f | %-15.4f | %-10.4f\n', ...
            SNR_dB(i_snr), PER_awgn(i_snr), PER_rayl(i_snr), PER_tdl(i_snr));
end

%% ── Построение графика ───────────────────────────────────────────────────
fig = figure('Color', 'white', ...
             'Units', 'centimeters', ...
             'Position', [3 3 16 12]);

ax = axes('Parent', fig, ...
          'Color',     'white', ...
          'FontName',  'Times New Roman', ...
          'FontSize',  12, ...
          'Box',       'on', ...
          'YScale',    'log', ...
          'XGrid',     'on', ...
          'YGrid',     'on', ...
          'GridColor', [0.85 0.85 0.85], ...
          'GridAlpha', 1, ...
          'TickDir',   'in', ...
          'LineWidth', 0.8);
hold(ax, 'on');

mstep = 4;
idx_m = 1 : mstep : length(SNR_dB);

% АБГШ
h1 = plot(ax, SNR_dB, max(PER_awgn, 1e-3), ...
          'k-', 'LineWidth', 1.5);
plot(ax, SNR_dB(idx_m), max(PER_awgn(idx_m), 1e-3), ...
     'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineStyle', 'none');

% Релей плоский
c_blue = [0.13 0.40 0.67];
h2 = plot(ax, SNR_dB, max(PER_rayl, 1e-3), ...
          '--', 'Color', c_blue, 'LineWidth', 1.5);
plot(ax, SNR_dB(idx_m), max(PER_rayl(idx_m), 1e-3), ...
     's', 'Color', c_blue, 'MarkerFaceColor', c_blue, ...
     'MarkerSize', 5, 'LineStyle', 'none');

% Релей TDL
c_red = [0.70 0.09 0.09];
h3 = plot(ax, SNR_dB, max(PER_tdl, 1e-3), ...
          ':', 'Color', c_red, 'LineWidth', 1.5);
plot(ax, SNR_dB(idx_m), max(PER_tdl(idx_m), 1e-3), ...
     '^', 'Color', c_red, 'MarkerFaceColor', c_red, ...
     'MarkerSize', 5, 'LineStyle', 'none');

% Горизонтальная линия diversity floor
yline(ax, 0.05, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
text(ax, -13, 0.065, 'Diversity floor ≈ 0.05', ...
     'FontName', 'Times New Roman', 'FontSize', 9, 'Color', [0.5 0.5 0.5]);

% Вертикальная линия порога чувствительности
xline(ax, SNR_thr_dB, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
text(ax, SNR_thr_dB + 0.2, 1.5e-3, ...
     sprintf('%.1f дБ', SNR_thr_dB), ...
     'FontName', 'Times New Roman', 'FontSize', 9, 'Color', [0.5 0.5 0.5]);

%% ── Подписи ─────────────────────────────────────────────────────────────
xlabel(ax, 'ОСШ (SNR), дБ', ...
       'FontName', 'Times New Roman', 'FontSize', 12);
ylabel(ax, 'Вероятность потери пакета (PER)', ...
       'FontName', 'Times New Roman', 'FontSize', 12);

legend(ax, [h1 h2 h3], {'АБГШ', 'Релей (плоский)', 'Релей TDL (4 луча)'}, ...
       'Location', 'southwest', ...
       'FontName', 'Times New Roman', ...
       'FontSize', 10, ...
       'Box',      'on');

xlim(ax, [min(SNR_dB) max(SNR_dB)]);
ylim(ax, [1e-3 1]);
xticks(ax, min(SNR_dB):2:max(SNR_dB));

%% ── Сохранение ───────────────────────────────────────────────────────────
set(fig, 'PaperPositionMode', 'auto');
print(fig, 'fig_PER_channel_comparison', '-dpng', '-r300');
fprintf('[OK] Сохранено: fig_PER_channel_comparison.png\n');

end

%% ── Вспомогательная функция: детектор ────────────────────────────────────
function sym = find_peak(spectrum)
% Возвращает индекс максимального бина (0-индексирован)
    [~, idx] = max(abs(spectrum));
    sym = idx - 1;   % переход к 0-индексации
end
