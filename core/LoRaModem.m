classdef LoRaModem < handle
    properties
        phy
        payloadLenBits

        % Эталонные настройки PHY
        cr_init
        has_header_init
        crc_init
        preamble_len_init
        rx_ignore_crc
        ref_chirp_
        cfg
    end

    methods
        function obj = LoRaModem(rf_freq, sf, bw, fs, varargin)
            p = inputParser;
            addParameter(p, 'CR', 4);
            addParameter(p, 'HasHeader', true);
            addParameter(p, 'UseCRC', true);
            addParameter(p, 'PreambleLen', 8);
            addParameter(p, 'FastMode', false);
            addParameter(p, 'RxIgnoreCRC', false);
            parse(p, varargin{:});

            obj.cr_init           = p.Results.CR;
            obj.has_header_init   = p.Results.HasHeader;
            obj.crc_init          = p.Results.UseCRC;
            obj.preamble_len_init = p.Results.PreambleLen;
            obj.rx_ignore_crc     = logical(p.Results.RxIgnoreCRC);

            obj.phy = LoRaPHY(rf_freq, sf, bw, fs);
            obj.cfg = makeLoRaConfig(rf_freq, sf, bw, fs);
            obj.phy.cr           = obj.cr_init;
            obj.phy.has_header   = obj.has_header_init;
            obj.phy.crc          = obj.crc_init;
            obj.phy.preamble_len = obj.preamble_len_init;
            obj.phy.fast_mode    = p.Results.FastMode;
            obj.phy.is_debug     = false;

            obj.payloadLenBits = [];
        end

        function [txSig, payloadBytes, nBitsUsed] = modulate(obj, bitsIn)
            if size(bitsIn, 2) > 1
                bitsIn = bitsIn(:);
            end

            obj.payloadLenBits = numel(bitsIn);

            % Восстанавливаем эталонные настройки PHY перед каждым TX
            obj.phy.cr           = obj.cr_init;
            obj.phy.has_header   = obj.has_header_init;
            obj.phy.crc          = obj.crc_init;
            obj.phy.preamble_len = obj.preamble_len_init;

            nBits = numel(bitsIn);
            nPad  = mod(8 - mod(nBits, 8), 8);
            if nPad == 8, nPad = 0; end
            bitsPadded = [bitsIn; zeros(nPad, 1, 'like', bitsIn)];
            nBitsUsed  = numel(bitsPadded);

            payloadBytes = obj.bits2bytes(bitsPadded);
            symbols      = obj.phy.encode(payloadBytes);
            txSig        = obj.phy.modulate(symbols);
        end

        function [bitsOut, rx_ok, crc_ok] = demodulate(obj, rxSig)
            % demodulate — демодулирует принятый сигнал.
            %
            % Выходные аргументы:
            %   bitsOut — восстановленные биты (может быть пустым при отказе)
            %   rx_ok   — true если декодер вернул данные нужной длины
            %   crc_ok  — true если CRC совпал (валидно только при rx_ok=true
            %             и UseCRC=true; при rx_ok=false всегда false)
            %
            % Матрица состояний:
            %   rx_ok=false, crc_ok=false → пакет полностью потерян (нет преамбулы /
            %                               невалидный заголовок / длина неверна)
            %   rx_ok=true,  crc_ok=false → данные получены, но CRC не сошёлся
            %                               (битовые ошибки есть, но бит сравнить можно)
            %   rx_ok=true,  crc_ok=true  → пакет принят успешно

            if size(rxSig, 2) > 1
                rxSig = rxSig(:);
            end

            bitsOut = false(0, 1);
            rx_ok   = false;
            crc_ok  = false;

            try
                [sym_rx, ~, ~] = obj.phy.demodulate(rxSig);

                if isempty(sym_rx)
                    return;   % нет преамбулы → полная потеря
                end

                % Временно снимаем CRC-проверку внутри decode(), чтобы
                % получить байты даже при несовпадении CRC. CRC проверим сами.
                crc_saved    = obj.phy.crc;
                obj.phy.crc  = false;   % decode без исключения при CRC-fail

                data_rx = [];
                try
                    [data_rx, ~] = obj.phy.decode(sym_rx);
                catch
                    obj.phy.crc = crc_saved;
                    return;   % невалидный заголовок или иная ошибка decode
                end

                obj.phy.crc = crc_saved;   % восстанавливаем

                if isempty(data_rx)
                    return;
                end

                % Конвертируем байты в биты и обрезаем до нужной длины
                bitsFull = obj.bytes2bits(uint8(data_rx));

                if ~isempty(obj.payloadLenBits) && numel(bitsFull) >= obj.payloadLenBits
                    bitsOut = bitsFull(1:obj.payloadLenBits);
                    rx_ok   = true;
                else
                    return;   % длина не та → потеря
                end

                % Проверка CRC вручную, если она включена в настройках.
                % LoRaPHY вычисляет CRC как последние 2 байта payload.
                % Сравниваем CRC переданных данных с CRC принятых данных.
                if crc_saved
                    crc_ok = obj.verifyCrc(data_rx, obj.phy);
                else
                    crc_ok = true;   % CRC отключён → считаем ок
                end

            catch
                bitsOut = false(0, 1);
                rx_ok   = false;
                crc_ok  = false;
            end
        end
        function [h_est, var_est, h_per_sym] = estimateChannel(obj, rxSig, snr_dB, s_ref)
        % estimateChannel  ML/MMSE-оценка коэффициента плоского замирающего канала
        %                  по пилотным символам преамбулы LoRa (Шаг 1 алгоритма
        %                  компенсации канальных искажений).
        %
        % ---------------------------------------------------------------
        % ФИЗИЧЕСКОЕ ОБОСНОВАНИЕ (Акимов/Бакут, гл. 1.2):
        %
        %   Модель наблюдения (observation model) для k-го символа преамбулы:
        %       r_k[n] = h · s_k[n] + w_k[n],   n = 0..Ns-1,  k = 1..N_p
        %
        %   где:
        %       h    — комплексный коэффициент плоского замирающего канала,
        %              h ~ CN(0, 1)  (замирания Релея)
        %       s_k  — известный k-й up-chirp преамбулы (пилотный символ)
        %       w_k  — АБГШ, w_k ~ CN(0, sigma_n^2 · I)
        %       Ns   — число отсчётов на символ: Ns = 2^SF · (fs/BW)
        %
        %   ML-оценка (maximum likelihood estimation) по одному символу:
        %       h_hat_k = <r_k, s_ref> / <s_ref, s_ref> = (r_k' · s_ref) / Ns
        %
        %   MMSE-усреднение по N_p символам преамбулы (оптимальное накопление):
        %       h_est = (1/N_p) · Σ_{k=1}^{N_p} h_hat_k
        %
        %   Теоретическая дисперсия (нижняя граница Крамера–Рао, CRLB):
        %       Var(h_est) = sigma_n^2 / (N_p · Ns)
        %                 = 1 / (N_p · Ns · SNR_lin)
        %
        %   Выигрыш накопления по N_p = 8 символам: -10·log10(8) ≈ -9 дБ к Var.
        %
        %   Примечание о нормировке SNR:
        %       SNR_lin задан как отношение мощности сигнала к мощности шума
        %       на один комплексный отсчёт (per-sample SNR convention).
        %       Эта нормировка соответствует SimpleChannel.addAwgn().
        % ---------------------------------------------------------------
        %
        % ВХОДНЫЕ АРГУМЕНТЫ:
        %   rxSig  [Mx1 complex] — принятый baseband-сигнал (column vector)
        %   snr_dB [scalar]      — ОСШ (SNR) в дБ (per-sample convention)
        %
        % ВЫХОДНЫЕ АРГУМЕНТЫ:
        %   h_est     [complex]  — MMSE-оценка коэффициента канала
        %   var_est   [double]   — теоретическая дисперсия: 1/(N_p·Ns·SNR_lin)
        %                          (используется как P_0 — начальная ковариация
        %                           для фильтра Калмана на Шаге 2)
        %   h_per_sym [N_p x 1]  — ML-оценки по отдельным символам преамбулы
        %                          (для анализа качества оценки и отладки)
        %
        % ПРИМЕР ИСПОЛЬЗОВАНИЯ (из сценария):
        %   [txSig, ~, ~]      = modem.modulate(bits_tx);
        %   rxSig              = channel.pass(txSig);
        %   [h_est, var_est]   = modem.estimateChannel(rxSig, snr_dB);
        %   % Инициализация фильтра Калмана (Шаг 2):
        %   h_kalman = h_est;   P_kalman = var_est;
        %
        % ЗАВИСИМОСТИ: нет внешних тулбоксов (только стандартные функции MATLAB)
        % ---------------------------------------------------------------
        
            %% --- Параметры дискретизации сигнала
            os = obj.cfg.os;    % oversampling factor (из LoRaConfig)
            N  = obj.cfg.N;     % chips per symbol   (из LoRaConfig)
            Ns = obj.cfg.Ns;    % samples per symbol (из LoRaConfig)

            %% --- Генерация референсного up-chirp (эталонного пилотного символа)
            %
            % Стандартный CSS up-chirp с линейным ЛЧМ от -BW/2 до +BW/2:
            %   f_inst[n] = -BW/2 + BW · n/Ns,   n = 0..Ns-1
            %   phi[n]    = (2π/fs) · Σ_{k=0}^{n} f_inst[k]
            %             = π · n · (n/Ns - 1/os)
            %
            if nargin < 4 || isempty(s_ref)
                s_ref = obj.getRefChirp();   % fallback: из своего экземпляра phy
            end
            s_energy = real(s_ref' * s_ref);
        
            %% --- Определение числа используемых символов преамбулы
            N_p     = obj.preamble_len_init;                       % длина преамбулы (по умолчанию 8)
            N_avail = floor(length(rxSig) / Ns);            % доступных символов в rxSig
            N_use   = min(N_p, N_avail);
        
            % Граничный случай: сигнал слишком короткий
            if N_use < 1
                warning('LoRaModem:estimateChannel:tooShort', ...
                    ['estimateChannel: rxSig слишком короткий. ' ...
                     'Ожидалось не менее %d отсчётов (1 символ), получено %d.'], ...
                    Ns, length(rxSig));
                h_est     = 1 + 0j;        % нейтральная оценка
                var_est   = Inf;           % бесконечная дисперсия — оценка ненадёжна
                h_per_sym = h_est;
                return;
            end
        
            if N_use < N_p
                warning('LoRaModem:estimateChannel:partialPreamble', ...
                    ['estimateChannel: используется %d из %d символов преамбулы ' ...
                     '(rxSig коротковат для полной преамбулы).'], N_use, N_p);
            end
        
            %% --- ML-оценка по каждому символу преамбулы
            h_per_sym = complex(zeros(N_use, 1));
        
            for k = 1 : N_use
                % Вырезаем k-й символ из принятого сигнала
                idx       = (k-1)*Ns + 1 : k*Ns;
                r_k       = rxSig(idx);                     % принятый k-й символ
        
                % Корреляция с референсным chirp: <r_k, s_ref> / <s_ref, s_ref>
                h_per_sym(k) = (s_ref' * r_k) / s_energy;
            end
        
            %% --- MMSE-усреднение по всем использованным символам
            %
            % При независимом аддитивном шуме w_k оценки h_hat_k некоррелированы,
            % поэтому усреднение является оптимальным (minimum variance estimator):
            %   h_est = (1/N_use) · Σ h_hat_k
            h_est = mean(h_per_sym);
        
            %% --- Теоретическая дисперсия оценки (CRLB)
            %
            % sigma_n^2 = 1/SNR_lin  (per-sample noise power)
            % Var(h_est) = sigma_n^2 / (N_use · Ns) = 1 / (N_use · Ns · SNR_lin)
            snr_lin = 10^(snr_dB / 10);
            var_est = 1 / (N_use * Ns * snr_lin);
        
        end
        function s_ref = getRefChirp(obj)
        % getRefChirp  Возвращает эталонный up-chirp символ, извлечённый из
        %              реального TX-сигнала LoRaPHY (не аналитическую аппроксимацию).
        %
        % При первом вызове генерирует минимальный LoRa-пакет, извлекает первый
        % символ преамбулы и сохраняет в кэше (obj.ref_chirp_).
        % При повторных вызовах возвращает кэшированное значение — без лишних
        % вызовов phy.encode() / phy.modulate() в цикле Монте-Карло.
        
            if isempty(obj.ref_chirp_)
                Ns = obj.cfg.Ns;    % samples per symbol (из LoRaConfig)

                % Сохраняем payloadLenBits — modulate() его перезапишет
                saved_len = obj.payloadLenBits;
        
                % Генерируем минимальный пакет (8 байт = 64 бита)
                % Содержимое payload не важно — нужен только первый символ преамбулы
                [txRef, ~, ~] = obj.modulate(false(64, 1));
        
                % Восстанавливаем состояние
                obj.payloadLenBits = saved_len;
        
                % Первый символ преамбулы — это и есть эталонный up-chirp
                obj.ref_chirp_ = txRef(1 : Ns);
            end
        
            s_ref = obj.ref_chirp_;
        end
    end

    methods (Access = private)
        function ok = verifyCrc(obj, data_rx, phy)
            % verifyCrc — сравнивает CRC принятого payload с вычисленным.
            %
            % LoRaPHY кладёт CRC как два последних байта в data_rx
            % (при has_header=true и crc=true). Вычисляем CRC по
            % первым (end-2) байтам и сравниваем с последними двумя.
            %
            % Если данных недостаточно для CRC-проверки → считаем fail.
            ok = false;
            try
                n = numel(data_rx);
                if n < 3
                    return;
                end
                payload      = data_rx(1:end-2);
                crc_received = data_rx(end-1:end);
                crc_calc     = phy.calc_crc(payload);
                ok = isequal(uint8(crc_received), uint8(crc_calc));
            catch
                ok = false;
            end
        end
    end

    methods (Static)
        function bytes = bits2bytes(bits)
            nBits = numel(bits);
            if mod(nBits, 8) ~= 0
                error('bits2bytes: число бит должно быть кратно 8');
            end
            bits        = reshape(bits, 8, []).';
            bits        = double(bits);
            powers      = 2.^(0:7);
            bytesDouble = bits * powers.';
            bytes       = uint8(bytesDouble);
        end

        function bits = bytes2bits(bytes)
            bytes  = uint8(bytes(:));
            nBytes = numel(bytes);
            bits   = false(nBytes * 8, 1);
            for k = 1:nBytes
                val = bytes(k);
                for b = 0:7
                    bits((k-1)*8 + b + 1) = bitget(val, b+1);
                end
            end
        end
    end
end
