classdef LoRaSimulator < handle
    properties
        modem      % объект LoRaModem
        channel    % объект SimpleChannel
    end

    methods
        function obj = LoRaSimulator(modem, channel)
            % Конструктор
            obj.modem   = modem;
            obj.channel = channel;
        end

        function [ber, per] = run(obj, Npkts, payloadLenBits)
            % Метрики:
            % per = доля erasure (payload не восстановлен до нужной длины)
            % ber = BER по payload среди успешно восстановленных кадров (rx_ok == true)
        
            nBitErr   = 0;
            nPktErr   = 0;   % erasures
            totalBits = 0;
        
            for k = 1:Npkts
                bits_tx = randi([0 1], payloadLenBits, 1, 'logical');
        
                [txSig, ~, ~] = obj.modem.modulate(bits_tx);
        
                rxSig = obj.channel.pass(txSig);
        
                [bits_rx, rx_ok] = obj.modem.demodulate(rxSig);
        
                if rx_ok
                    % payload должен быть ровно нужной длины
                    if numel(bits_rx) ~= payloadLenBits
                        nPktErr = nPktErr + 1;
                        continue;
                    end
        
                    nBitErr   = nBitErr + sum(bits_tx ~= bits_rx);
                    totalBits = totalBits + payloadLenBits;
                else
                    nPktErr = nPktErr + 1;
                end
            end
        
            if totalBits == 0
                ber = NaN;   % не было ни одного принятого payload
            else
                ber = nBitErr / totalBits;
            end
        
            per = nPktErr / Npkts;
        end
        % function [ber, per] = run(obj, Npkts, payloadLenBits)
        %     % Запускает имитацию пакетов указанной длины (в битах)
        %     nBitErr   = 0;
        %     nPktErr   = 0;
        %     totalBits = 0;
        % 
        %     for k = 1:Npkts
        %         % Генерируем случайные биты
        %         bits_tx = randi([0 1], payloadLenBits, 1, 'logical');
        % 
        %         % Модуляция
        %         [txSig, ~, ~] = obj.modem.modulate(bits_tx);
        % 
        %         % Канал
        %         rxSig = obj.channel.pass(txSig);
        % 
        %         % Демодуляция
        %         [bits_rx, ~] = obj.modem.demodulate(rxSig);
        % 
        %         % Подсчет битовых ошибок
        %         L = min(numel(bits_tx), numel(bits_rx));
        %         if L > 0
        %             nBitErr   = nBitErr + sum(bits_tx(1:L) ~= bits_rx(1:L));
        %             totalBits = totalBits + L;
        %         end
        % 
        %         % Пакет считаем успешным, если RX не короче TX и все биты совпали
        %         if numel(bits_rx) >= numel(bits_tx)
        %             pkt_ok = all(bits_tx == bits_rx(1:numel(bits_tx)));
        %         else
        %             pkt_ok = false;
        %         end
        % 
        %         % % Отладочный вывод для первого пакета
        %         % if k == 1
        %         %     disp('TX bits:'); disp(bits_tx.');
        %         %     disp('RX bits:'); disp(bits_rx.');
        %         %     fprintf('pkt_ok = %d\n', pkt_ok);
        %         % end
        % 
        %         if ~pkt_ok
        %             nPktErr = nPktErr + 1;
        %         end
        %     end
        % 
        %     % 6) Финальные метрики
        %     if totalBits == 0
        %         ber = NaN;   % на случай, если вообще ничего не сравнили
        %     else
        %         ber = nBitErr / totalBits;
        %     end
        % 
        %     per = nPktErr / Npkts;
        % end
    end
end
