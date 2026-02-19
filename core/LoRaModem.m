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
