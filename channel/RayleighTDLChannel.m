classdef RayleighTDLChannel < SimpleChannel
    % RayleighTDLChannel — частотно-селективный канал Релея, TDL-модель.
    %
    % Поддерживает два режима замираний:
    %
    %   'PerPacket'  — один коэффициент h на весь пакет (block fading).
    %                  Устанавливается автоматически, если SF/BW не заданы.
    %                  Поведение идентично предыдущей версии класса.
    %
    %   'PerSymbol'  — коэффициент h обновляется раз в символ по модели AR(1).
    %                  Устанавливается автоматически при задании SF и BW.
    %                  Физически корректно: T_c > T_sym (медленные замирания),
    %                  канал постоянен внутри символа, но меняется между ними.
    %
    % Модель AR(1) для каждого луча k:
    %
    %   h_k[n] = rho * h_k[n-1] + sqrt(1 - rho^2) * sigma_k * w[n]
    %   w[n]   ~ CN(0, 1)
    %   rho    = J0(2*pi * f_D * T_sym)
    %
    % При rho=1 (f_D=0) канал замерзает — эквивалент PerPacket на масштабе пакета.
    % При rho=0         — независимые символы (AR(0)).
    % При rho≈0.92      — реалистичный сценарий: v=30 м/с, SF7, f_c=868 МГц.
    %
    % Корреляция AR(1) воспроизводит корреляционную функцию модели Джейкса
    % на временно́м масштабе символа без использования Communications Toolbox.
    %
    % Обратная совместимость:
    %   Все существующие вызовы вида
    %       RayleighTDLChannel(fs, snr_dB, cfo_Hz, 'PathDelays', ..., 'Seed', ...)
    %   работают без изменений — режим PerPacket включается автоматически.
    %
    % Использование (PerSymbol, рекомендуется):
    %   ch = RayleighTDLChannel(fs, snr_dB, cfo_Hz, ...
    %       'SF', 7, 'BW', 125e3, ...
    %       'DopplerHz', 87, ...
    %       'PathDelays', [0, 0.5e-6, 1.5e-6], ...
    %       'PathGains',  [0, -6, -12], ...
    %       'Seed', 42);
    %
    % Использование (PerPacket, для сравнения):
    %   ch = RayleighTDLChannel(fs, snr_dB, cfo_Hz, ...
    %       'PathDelays', [0, 0.5e-6, 1.5e-6], ...
    %       'PathGains',  [0, -6, -12], ...
    %       'Seed', 42);
    %
    % Совместимость: MATLAB R2024a, без дополнительных тулбоксов.

    % ------------------------------------------------------------------
    % Публичные свойства
    % ------------------------------------------------------------------
    properties
        pathDelays_s        % [L×1] задержки лучей, с
        pathGains_dB        % [L×1] средние мощности лучей, дБ
        seed                % seed RNG ([] = случайный)
        sf                  % Spreading Factor ([] если не задан)
        bw                  % полоса сигнала, Гц ([] если не задана)
        dopplerHz           % максимальная допплеровская частота, Гц
        fadingGranularity   % 'PerPacket' | 'PerSymbol'
    end

    % ------------------------------------------------------------------
    % Приватные свойства — внутреннее состояние
    % ------------------------------------------------------------------
    properties (Access = private)
        h_current       % [L×1] текущее состояние AR(1) (инициализируется в конструкторе)
        h_initial       % [L×1] начальное состояние для resetRng()
        N_sym           % длина символа в отсчётах = round(2^sf / bw * fs)
        rho             % коэффициент AR(1): J0(2*pi * f_D * T_sym)
        gains_lin_norm  % [L×1] нормированные линейные мощности лучей (предвычислено)
        rng_state       % состояние RNG при инициализации
    end

    % ------------------------------------------------------------------
    % Публичные методы
    % ------------------------------------------------------------------
    methods

        % --------------------------------------------------------------
        % Конструктор
        % --------------------------------------------------------------
        function obj = RayleighTDLChannel(fs, snr_dB, cfo_Hz, varargin)
            % Вызов конструктора базового класса
            obj@SimpleChannel(fs, snr_dB, cfo_Hz);

            % Параметры TDL по умолчанию (3 луча, tau_rms ≈ 0.6 мкс)
            p = inputParser;
            addParameter(p, 'PathDelays', [0, 0.5e-6, 1.5e-6]);
            addParameter(p, 'PathGains',  [0, -6,     -12]);
            addParameter(p, 'Seed',       []);
            addParameter(p, 'SF',         []);   % задать для режима PerSymbol
            addParameter(p, 'BW',         []);   % задать для режима PerSymbol
            addParameter(p, 'DopplerHz',  0);    % f_D, Гц
            addParameter(p, 'FadingGranularity', []);  % [] = автоопределение
            parse(p, varargin{:});

            obj.pathDelays_s  = p.Results.PathDelays(:);
            obj.pathGains_dB  = p.Results.PathGains(:);
            obj.seed          = p.Results.Seed;
            obj.sf            = p.Results.SF;
            obj.bw            = p.Results.BW;
            obj.dopplerHz     = p.Results.DopplerHz;

            % Проверка согласованности профиля
            assert(numel(obj.pathDelays_s) == numel(obj.pathGains_dB), ...
                'RayleighTDLChannel: PathDelays и PathGains должны иметь одинаковую длину.');

            % Автоматическое определение режима замираний:
            %   SF+BW заданы → PerSymbol (физически корректно)
            %   иначе        → PerPacket (обратная совместимость)
            fg = p.Results.FadingGranularity;
            if isempty(fg)
                if ~isempty(obj.sf) && ~isempty(obj.bw)
                    fg = 'PerSymbol';
                else
                    fg = 'PerPacket';
                end
            end
            obj.fadingGranularity = fg;

            % Проверка: PerSymbol требует SF и BW
            if strcmp(obj.fadingGranularity, 'PerSymbol')
                assert(~isempty(obj.sf) && ~isempty(obj.bw), ...
                    ['RayleighTDLChannel: для FadingGranularity=''PerSymbol'' ', ...
                     'необходимо задать параметры ''SF'' и ''BW''.']);
            end

            % Вычисление параметров символа и AR(1)
            if ~isempty(obj.sf) && ~isempty(obj.bw)
                obj.N_sym = round((2^obj.sf / obj.bw) * obj.fs);
                T_sym     = 2^obj.sf / obj.bw;
                % Коэффициент AR(1): rho = J0(2*pi * f_D * T_sym)
                % J0 — функция Бесселя первого рода нулевого порядка (besselj)
                obj.rho   = besselj(0, 2*pi * obj.dopplerHz * T_sym);
            else
                obj.N_sym = [];
                obj.rho   = 0;
            end

            % Предвычисление нормированных линейных мощностей лучей.
            % Нормировка: sum(gains_lin_norm) = 1 → канал не вносит
            % среднего усиления/ослабления.
            gains        = 10.^(obj.pathGains_dB / 10);
            obj.gains_lin_norm = gains / sum(gains);

            % Инициализация генератора случайных чисел
            if ~isempty(obj.seed)
                rng(obj.seed, 'twister');
            end
            obj.rng_state = rng;

            % Инициализация состояния AR(1) из стационарного распределения.
            % Это физически корректно: h_k[0] ~ CN(0, gains_norm(k)).
            obj.h_current = obj.genRayleighCoeffs();
            obj.h_initial = obj.h_current;
        end

        % --------------------------------------------------------------
        % Прохождение сигнала через канал
        % --------------------------------------------------------------
        function y = pass(obj, x)
            % pass — применяет TDL-замирания, CFO и АБГШ.
            %
            % Цепочка: applyTDL → CFO-поворот → AWGN
            % Длина выхода равна длине входа (критично для детектора преамбулы).

            if size(x, 2) > 1
                x = x(:);
            end
            P_tx = mean(abs(x).^2);

            % PerPacket: обновить h независимо перед каждым пакетом.
            % PerSymbol: h_current — это AR(1)-состояние, оно обновляется
            %            внутри applyTDL(), обеспечивая непрерывность процесса
            %            между пакетами.
            if strcmp(obj.fadingGranularity, 'PerPacket')
                obj.h_current = obj.genRayleighCoeffs();
            end

            if strcmp(obj.fadingGranularity, 'PerPacket') || obj.rho >= 1.0
                % При rho=1 (v=0) PerSymbol вырождается в канал полностью в глубоком замирании.
                % Регенерируем h независимо для каждого пакета, как в PerPacket.
                obj.h_current = obj.genRayleighCoeffs();
            end

            % Применение TDL (замирания + многолучевость)
            x_mp = obj.applyTDL(x);

            % CFO-поворот
            t     = (0:length(x_mp)-1).' / obj.fs;
            x_cfo = x_mp .* exp(1j * 2 * pi * obj.cfo_Hz .* t);

            % АБГШ (SNR нормируется по мощности переданного сигнала)
            snr_lin = 10^(obj.snr_dB / 10);
            P_noise = P_tx / snr_lin;
            noise   = sqrt(P_noise / 2) * (randn(size(x_cfo)) + 1j * randn(size(x_cfo)));
            y       = x_cfo + noise;
        end

        % --------------------------------------------------------------
        % Вычисление tau_rms для текущего профиля
        % --------------------------------------------------------------
        function tau = calcTauRms(obj)
            % calcTauRms — СКО задержек по профилю PDP, секунды.
            tau_mean = sum(obj.pathDelays_s .* obj.gains_lin_norm);
            tau      = sqrt(sum((obj.pathDelays_s - tau_mean).^2 .* obj.gains_lin_norm));
        end

        % --------------------------------------------------------------
        % Вывод параметров канала
        % --------------------------------------------------------------
        function printInfo(obj)
            fprintf('--- RayleighTDLChannel ---\n');
            fprintf('  SNR        = %.1f дБ\n', obj.snr_dB);
            fprintf('  CFO        = %.1f Гц\n', obj.cfo_Hz);
            fprintf('  fs         = %.3g Гц\n', obj.fs);
            fprintf('  Режим замираний : %s\n', obj.fadingGranularity);
            fprintf('  Число лучей: %d\n', numel(obj.pathDelays_s));
            fprintf('  %-6s %-16s %-14s\n', 'Луч', 'Задержка, мкс', 'Мощность, дБ');
            for k = 1:numel(obj.pathDelays_s)
                fprintf('  %-6d %-16.3f %-14.1f\n', k, ...
                    obj.pathDelays_s(k)*1e6, obj.pathGains_dB(k));
            end
            tau_rms = obj.calcTauRms();
            fprintf('  tau_rms    = %.3f мкс\n', tau_rms*1e6);
            fprintf('  Bc (1/5·tau_rms) = %.1f кГц\n', ...
                1 / (5 * max(tau_rms, eps)) / 1e3);
            if ~isempty(obj.sf)
                T_sym = 2^obj.sf / obj.bw;
                fprintf('  SF=%d, BW=%.0f кГц\n', obj.sf, obj.bw/1e3);
                fprintf('  T_sym = %.3f мс,  N_sym = %d отсчётов\n', ...
                    T_sym*1e3, obj.N_sym);
                fprintf('  f_D   = %.1f Гц\n', obj.dopplerHz);
                fprintf('  rho   = %.4f  (J0(2π·f_D·T_sym))\n', obj.rho);
                rho_pkt = besselj(0, 2*pi * obj.dopplerHz * ...
                    (obj.N_sym / obj.fs) * ceil(2000 / obj.N_sym));
                fprintf('  Замечание: T_c/T_sym ≈ %.1f → %s\n', ...
                    0.423 / (max(obj.dopplerHz, 0.1) * T_sym), ...
                    ternary(strcmp(obj.fadingGranularity, 'PerSymbol'), ...
                        'символьные замирания (физически корректно)', ...
                        'блоковые замирания'));
            end
            if ~isempty(obj.seed)
                fprintf('  Seed       = %d\n', obj.seed);
            else
                fprintf('  Seed       = случайный\n');
            end
            fprintf('--------------------------\n');
        end

        % --------------------------------------------------------------
        % Сброс генератора к начальному состоянию
        % --------------------------------------------------------------
        function resetRng(obj)
            % resetRng — восстанавливает состояние RNG и h_current
            % до значений на момент конструктора. Позволяет повторить
            % ровно тот же канал при повторном запуске симуляции.
            rng(obj.rng_state);
            obj.h_current = obj.h_initial;
        end

    end % methods (public)

    % ------------------------------------------------------------------
    % Приватные методы
    % ------------------------------------------------------------------
    methods (Access = private)

        % --------------------------------------------------------------
        % Генерация вектора комплексных коэффициентов Релея
        % --------------------------------------------------------------
        function h = genRayleighCoeffs(obj)
            % genRayleighCoeffs — одна реализация h из стационарного
            % распределения: h_k ~ CN(0, gains_norm(k)).
            %
            % Используется:
            %   - в конструкторе (инициализация состояния AR(1));
            %   - в pass() при режиме PerPacket (независимый h на пакет).
            L = numel(obj.pathDelays_s);
            h = zeros(L, 1);
            for k = 1:L
                sigma = sqrt(obj.gains_lin_norm(k) / 2);
                h(k)  = sigma * (randn + 1j * randn);
            end
        end

        % --------------------------------------------------------------
        % TDL-свёртка с замираниями
        % --------------------------------------------------------------
        function y = applyTDL(obj, x)
            % applyTDL — применяет TDL-канал к сигналу x.
            %
            % Дискретная модель (PerPacket):
            %   y[n] = sum_k h_k * x[n - d_k]
            %
            % Дискретная модель (PerSymbol):
            %   y[n] = sum_k h_k[sym(n)] * x[n - d_k]
            %
            % где sym(n) = floor(n / N_sym) — номер символа отсчёта n.
            %
            % Граничное условие: x[n] = 0 при n < 0.
            % Длина выхода = длине входа (хвост свёртки отбрасывается,
            % что при tau_max << T_sym несущественно).

            N = length(x);
            y = zeros(N, 1);

            if strcmp(obj.fadingGranularity, 'PerPacket')
                % ======================================================
                %  PerPacket: один h на весь пакет
                % ======================================================
                for k = 1:numel(obj.pathDelays_s)
                    d = round(obj.pathDelays_s(k) * obj.fs);
                    if d == 0
                        y = y + obj.h_current(k) * x;
                    elseif d < N
                        y(d+1:end) = y(d+1:end) + obj.h_current(k) * x(1:end-d);
                    end
                    % d >= N → луч за окном пакета, игнорируем
                end

            else
                % ======================================================
                %  PerSymbol: AR(1), один h на символ
                % ======================================================
                N_sym  = obj.N_sym;
                nSyms  = ceil(N / N_sym);
                L      = numel(obj.pathDelays_s);

                % Генерация последовательности коэффициентов h[n]
                % для всех символов пакета по рекуррентной формуле AR(1):
                %
                %   h_k[n] = rho * h_k[n-1] + sqrt(1-rho^2) * sigma_k * w_k[n]
                %   w_k[n] ~ CN(0, 1)
                %
                % Начальное состояние h_current сохраняется между пакетами,
                % обеспечивая непрерывность AR(1)-процесса.
                h_seq  = zeros(L, nSyms);
                h_prev = obj.h_current;
                coeff  = sqrt(1 - obj.rho^2);

                for n = 1:nSyms
                    h_new = zeros(L, 1);
                    for k = 1:L
                        sigma   = sqrt(obj.gains_lin_norm(k) / 2);
                        w       = sigma * (randn + 1j * randn);
                        h_new(k) = obj.rho * h_prev(k) + coeff * w;
                    end
                    h_seq(:, n) = h_new;
                    h_prev      = h_new;
                end

                % Обновляем состояние для следующего пакета
                obj.h_current = h_prev;

                % Применение TDL: для каждого луча — посимвольная свёртка.
                % Каждый отсчёт y[m] использует h_k того символа, которому
                % принадлежит m: sym = ceil(m / N_sym).
                for k = 1:L
                    d = round(obj.pathDelays_s(k) * obj.fs);
                    if d >= N
                        continue;
                    end

                    for n = 1:nSyms
                        % Диапазон выходных отсчётов символа n
                        dst_s = (n-1)*N_sym + 1;
                        dst_e = min(n*N_sym, N);
                        h_k   = h_seq(k, n);

                        % Исходные отсчёты (с учётом задержки d)
                        src_s = dst_s - d;
                        src_e = dst_e - d;

                        if src_e < 1
                            continue;   % весь диапазон до начала сигнала
                        end

                        % Обрезаем левый край, если src_s < 1
                        skip  = max(0, 1 - src_s);
                        dst_s = dst_s + skip;
                        src_s = src_s + skip;

                        if dst_s > dst_e
                            continue;
                        end

                        y(dst_s:dst_e) = y(dst_s:dst_e) + h_k * x(src_s:src_e);
                    end
                end
            end
        end

    end % methods (private)

end

% Вспомогательная функция: тернарный оператор для printInfo()
function r = ternary(cond, a, b)
    if cond, r = a; else, r = b; end
end