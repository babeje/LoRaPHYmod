classdef LoRaSimulator < handle
    properties
        modem   % объект LoRaModem
        channel % объект SimpleChannel
    end
    
    methods
        function obj = LoRaSimulator(modem, channel)
            obj.modem   = modem;
            obj.channel = channel;
        end
        
        function [ber, per] = run(obj, Npkts, payloadLenBits)
            nBitErr   = 0;
            nPktErr   = 0;
            totalBits = 0;
            
            for k = 1:Npkts
                % Генерация случайных бит
                bits_tx = randi([0 1], payloadLenBits, 1, 'logical');
                
                % Модуляция
                [txSig, ~, ~] = obj.modem.modulate(bits_tx);
                
                % Проход через канал
                rxSig = obj.channel.pass(txSig);
                
                % Демодуляция
                [bits_rx, ok] = obj.modem.demodulate(rxSig);
                
                % Подсчет ошибок
                L = min(numel(bits_tx), numel(bits_rx));
                nBitErr   = nBitErr + sum(bits_tx(1:L) ~= bits_rx(1:L));
                totalBits = totalBits + L;
                
                if (~ok) || any(bits_tx(1:L) ~= bits_rx(1:L))
                    nPktErr = nPktErr + 1;
                end
            end
            
            ber = nBitErr / totalBits;
            per = nPktErr / Npkts;
        end
    end
end
