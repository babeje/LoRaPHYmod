function fig_BER_SF_sweep()
% fig_BER_SF_sweep  Зависимость BER от ОСШ для SF=7..12 в канале АБГШ
%
% Аналитическая оценка на основе границы объединения (union bound)
% для некогерентного CSS-детектора LoRa:
%
%   P_s ≈ (M-1)/2 · exp(−M·γ_BW/2),   M = 2^SF
%   BER = P_s / SF                      (равномерное распределение ошибок)
%   γ_BW = 10^(SNR_dB/10)              (ОСШ в полосе BW)
%
% Использование:
%   fig_BER_SF_sweep()
%
% Выход:
%   fig_BER_SF_sweep.png  (300 DPI, академический стиль)

%% ── Параметры ─────────────────────────────────────────────────────────────
SF_vals = 7:12;
SNR_dB  = (-25 : 0.25 : 5)';          % диапазон ОСШ, дБ

%% ── Аналитический расчёт BER ──────────────────────────────────────────────
SNR_lin = 10 .^ (SNR_dB / 10);
BER_mat = zeros(length(SNR_dB), length(SF_vals));

for k = 1:length(SF_vals)
    SF  = SF_vals(k);
    M   = 2 ^ SF;
    Ps  = (M - 1) / 2 .* exp(-M .* SNR_lin ./ 2);   % граница объединения
    BER = Ps ./ SF;
    BER = min(BER, 0.5);     % физический предел (случайное решение)
    BER = max(BER, 1e-7);    % нижний предел для отображения
    BER_mat(:, k) = BER;
end

%% ── Настройки оформления ──────────────────────────────────────────────────
colors  = { [0   0   0  ],   ...   % SF=7  чёрный
            [0.13 0.40 0.67], ...   % SF=8  тёмно-синий
            [0.70 0.09 0.09], ...   % SF=9  тёмно-красный
            [0.30 0.62 0.15], ...   % SF=10 тёмно-зелёный
            [0.49 0.25 0.0 ], ...   % SF=11 коричневый
            [0.37 0.24 0.60] };     % SF=12 тёмно-фиолетовый

lstyles = {'-', '--', ':', '-.', '--', ':'};
markers = {'o', 's',  '^', 'd',  'v',  'p'};
mstep   = 10;                        % каждые 10 точек — маркер

%% ── Построение ───────────────────────────────────────────────────────────
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

h_lines = gobjects(length(SF_vals), 1);

for k = 1:length(SF_vals)
    h_lines(k) = plot(ax, SNR_dB, BER_mat(:, k), ...
        'Color',     colors{k}, ...
        'LineStyle', lstyles{k}, ...
        'LineWidth', 1.5);

    idx_m = 1 : mstep : length(SNR_dB);
    plot(ax, SNR_dB(idx_m), BER_mat(idx_m, k), ...
        markers{k}, ...
        'Color',           colors{k}, ...
        'MarkerFaceColor', colors{k}, ...
        'MarkerSize',      5, ...
        'LineStyle',       'none');
end

%% ── Подписи осей и легенда ────────────────────────────────────────────────
xlabel(ax, 'ОСШ (SNR), дБ', ...
       'FontName', 'Times New Roman', 'FontSize', 12);
ylabel(ax, 'Вероятность битовой ошибки (BER)', ...
       'FontName', 'Times New Roman', 'FontSize', 12);

leg_str = arrayfun(@(sf) sprintf('SF = %d', sf), SF_vals, ...
                   'UniformOutput', false);
legend(ax, h_lines, leg_str, ...
       'Location', 'southwest', ...
       'FontName', 'Times New Roman', ...
       'FontSize', 10, ...
       'Box',      'on');

xlim(ax, [-25  5]);
ylim(ax, [1e-4 1]);
xticks(ax, -25:5:5);

%% ── Сохранение ───────────────────────────────────────────────────────────
set(fig, 'PaperPositionMode', 'auto');
print(fig, 'fig_BER_SF_sweep', '-dpng', '-r300');
fprintf('[OK] Сохранено: fig_BER_SF_sweep.png\n');
end
