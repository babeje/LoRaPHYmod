classdef LoRaModem < handle
    properties
        phy            % объект LoRaPHY
        payloadLenBits % длина полезной нагрузки в битах (для проверки)
    end
    
    methods
        function obj = LoRaModem(rf_freq, sf, bw, fs, varargin)
            % Параметры по умолчанию
            p = inputParser;
            addParameter(p, 'CR', 4);          % 1=4/5,2=4/6,3=4/7,4=4/8
            addParameter(p, 'HasHeader', true);
            addParameter(p, 'UseCRC', true);
            addParameter(p, 'PreambleLen', 8);
            addParameter(p, 'FastMode', false);
            parse(p, varargin{:});
            
            % Инициализация LoRaPHY
            obj.phy = LoRaPHY(rf_freq, sf, bw, fs);
            obj.phy.cr           = p.Results.CR;
            obj.phy.has_header   = p.Results.HasHeader;
            obj.phy.crc          = p.Results.UseCRC;
            obj.phy.preamble_len = p.Results.PreambleLen;
            obj.phy.fast_mode    = p.Results.FastMode;
            obj.phy.is_debug     = false;
            
            obj.payloadLenBits = [];   % можно задать позже
        end
        
        function [txSig, payloadBytes, nBitsUsed] = modulate(obj, bitsIn)
            % bitsIn: столбец 0/1
            if size(bitsIn,2) > 1
                bitsIn = bitsIn(:); % в столбец
            end
            
            obj.payloadLenBits = numel(bitsIn);
            
            % Дополняем до целого числа байт нулями
            nBits = numel(bitsIn);
            nPad  = mod(8 - mod(nBits,8), 8);
            if nPad == 8
                nPad = 0;
            end
            bitsPadded = [bitsIn; zeros(nPad,1,'like',bitsIn)];
            nBitsUsed  = numel(bitsPadded);
            
            % Биты -> байты
            payloadBytes = obj.bits2bytes(bitsPadded);
            
            % LoRaPHY TX
            symbols = obj.phy.encode(payloadBytes);
            txSig   = obj.phy.modulate(symbols);
        end
        
        function [bitsOut, ok] = demodulate(obj, rxSig)
            % LoRaPHY RX
            [sym_rx, ~, ~]   = obj.phy.demodulate(rxSig);
            [data_rx, chkOK] = obj.phy.decode(sym_rx);  % data_rx - полезная нагрузка
            
            % Байты -> биты
            bitsFull = obj.bytes2bits(data_rx);
            
            % Обрезаем до исходной длины
            if isempty(obj.payloadLenBits)
                bitsOut = bitsFull;
            else
                L = min(obj.payloadLenBits, numel(bitsFull));
                bitsOut = bitsFull(1:L);
            end
    
            % chkOK может быть не скаляром -> "сжимаем" до одного булева
            okChecksum = all(chkOK(:) == 1);
    
            % Условие "пакет принят"
            okLength = (numel(bitsOut) == obj.payloadLenBits);
            ok       = okChecksum && okLength;
        end
    end
    
    methods (Static)
        function bytes = bits2bytes(bits)
            % bits: столбец 0/1, длина кратна 8
            nBits = numel(bits);
            if mod(nBits,8) ~= 0
                error('bits2bytes: число бит должно быть кратно 8');
            end
    
            bits = reshape(bits, 8, []).';   % [Nbytes x 8]
    
            % Переводим в double, чтобы нормально перемножать
            bits = double(bits);             % 0/1 в double
    
            % Веса бит (LSB first): b0 + 2*b1 + ... + 128*b7
            powers = 2.^(0:7);               % double[1x8]
    
            % Матричное умножение даёт вектор байт (в double)
            bytesDouble = bits * powers.';   % [Nbytes x 1] double
    
            % Приводим к uint8
            bytes = uint8(bytesDouble);
        end
    
        function bits = bytes2bits(bytes)
            % bytes: вектор uint8 → столбец бит (LSB first)
            bytes  = uint8(bytes(:));     % в столбец
            nBytes = numel(bytes);
    
            bits = false(nBytes*8,1);     % логический столбец
    
            for k = 1:nBytes
                val = bytes(k);
                % биты 0..7 → позиции 1..8 (LSB first)
                for b = 0:7
                    bits((k-1)*8 + b + 1) = bitget(val, b+1);
                end
            end
        end
    end

    
    % methods (Static)
    %     function bytes = bits2bytes(bits)
    %         % bits: столбец 0/1, длина кратна 8
    %         nBits = numel(bits);
    %         if mod(nBits,8) ~= 0
    %             error('bits2bytes: число бит должно быть кратно 8');
    %         end
    %         bits = reshape(bits, 8, []).';
    %         % ЛСБ первым: b0 + 2*b1 + 4*b2 + ... + 128*b7
    %         powers = uint8(2.^(0:7));
    %         bytes  = uint8(bits) * powers.';
    %     end
    % 
    %     function bits = bytes2bits(bytes)
    %         % bytes: uint8
    %         bytes = uint8(bytes(:));
    %         nBytes = numel(bytes);
    %         bits = false(nBytes*8,1);
    %         for k = 1:nBytes
    %             val = bytes(k);
    %             for b = 0:7
    %                 bits((k-1)*8 + b + 1) = bitget(val, b+1); % снова ЛСБ первым
    %             end
    %         end
    %     end
    % end
end
