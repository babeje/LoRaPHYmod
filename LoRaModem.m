classdef LoRaModem < handle
    properties
        phy
        payloadLenBits
    
        % Эталонные настройки PHY
        cr_init
        has_header_init
        crc_init
        preamble_len_init
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
            
            % Сохраняем эталонные настройки
            obj.cr_init           = p.Results.CR;
            obj.has_header_init   = p.Results.HasHeader;
            obj.crc_init          = p.Results.UseCRC;
            obj.preamble_len_init = p.Results.PreambleLen;
        
            % Инициализируем LoRaPHY
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
            if size(bitsIn,2) > 1
                bitsIn = bitsIn(:);
            end
        
            obj.payloadLenBits = numel(bitsIn);
        
            % перед каждым TX возвращаем PHY к эталонным настройкам 
            obj.phy.cr           = obj.cr_init;
            obj.phy.has_header   = obj.has_header_init;
            obj.phy.crc          = obj.crc_init;
            obj.phy.preamble_len = obj.preamble_len_init;
        
            nBits = numel(bitsIn);
            nPad  = mod(8 - mod(nBits,8), 8);
            if nPad == 8
                nPad = 0;
            end
            bitsPadded = [bitsIn; zeros(nPad,1,'like',bitsIn)];
            nBitsUsed  = numel(bitsPadded);
        
            payloadBytes = obj.bits2bytes(bitsPadded);
            symbols      = obj.phy.encode(payloadBytes);
            txSig        = obj.phy.modulate(symbols);
        end
        
        function [bitsOut, ok] = demodulate(obj, rxSig)
            % LoRaPHY RX
            [sym_rx, ~, ~]   = obj.phy.demodulate(rxSig);
            [data_rx, chkVal] = obj.phy.decode(sym_rx);  % data_rx - полезная нагрузка
            
            % Байты -> биты
            bitsFull = obj.bytes2bits(data_rx);
            
            % Обрезаем до исходной длины
            if isempty(obj.payloadLenBits)
                bitsOut = bitsFull;
            else
                L = min(obj.payloadLenBits, numel(bitsFull));
                bitsOut = bitsFull(1:L);
            end

            % проверка по длине
            ok = (numel(bitsOut) == obj.payloadLenBits);
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
            % bytes: вектор uint8 в столбец бит (LSB first)
            bytes  = uint8(bytes(:));     % в столбец
            nBytes = numel(bytes);
    
            bits = false(nBytes*8,1);     % логический столбец
    
            for k = 1:nBytes
                val = bytes(k);
                % биты 0..7 в позиции 1..8 (LSB first)
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
