classdef LoRaSimulator < handle
    properties
        modem      % объект LoRaModem
        channel    % объект канала (SimpleChannel или потомок)
    end

    methods
        function obj = LoRaSimulator(modem, channel)
            obj.modem   = modem;
            obj.channel = channel;
        end

        function [ber, per, stats] = run(obj, Npkts, payloadLenBits)
            % run — основной цикл симуляции.
            %
            % Метрики:
            %   ber   — BER по всем пакетам (включая потерянные)
            %   per   — PER: доля пакетов с CRC-fail или полной потерей
            %   stats — структура с детальной статистикой
            %
            % Классификация исходов пакета:
            %
            %   [rx_ok=0, crc_ok=0] → ERASURE
            %     Декодер не смог восстановить данные нужной длины.
            %     Причина: нет преамбулы, невалидный заголовок, пустой decode.
            %     → PER: +1 (потеря)
            %     → BER: +0.5 на каждый бит (максимальная неопределённость)
            %
            %   [rx_ok=1, crc_ok=0] → CRC FAIL
            %     Байты получены нужной длины, но CRC не совпал.
            %     Данные есть — биты сравниваем побитово.
            %     → PER: +1 (пакет считается потерянным, как в реальном LoRa)
            %     → BER: считаем реальные битовые ошибки
            %
            %   [rx_ok=1, crc_ok=1] → SUCCESS
            %     Пакет принят корректно.
            %     → PER: 0
            %     → BER: считаем реальные битовые ошибки (обычно 0)
            %
            % Итоговый BER включает все три случая, что даёт
            % физически корректную кривую BER(SNR) без разрывов.

            % Счётчики
            nBitErr_real   = 0;   % реальные битовые ошибки (rx_ok=true)
            nBitErr_erased = 0;   % условные ошибки по erasure (0.5 * bits)
            nPktErr        = 0;   % потерянные пакеты (erasure + CRC fail)
            nErasure       = 0;   % полные потери (rx_ok=false)
            nCrcFail       = 0;   % приняты, но CRC не сошёлся
            nSuccess       = 0;   % успешно принятые пакеты

            for k = 1:Npkts
                bits_tx = logical(randi([0 1], payloadLenBits, 1));

                [txSig, ~, ~]        = obj.modem.modulate(bits_tx);
                rxSig                = obj.channel.pass(txSig);
                [bits_rx, rx_ok, crc_ok] = obj.modem.demodulate(rxSig);

                if ~rx_ok
                    % --- ERASURE: декодер не вернул данные нужной длины ---
                    % BER: не знаем ни одного бита → максимальная неопределённость
                    nBitErr_erased = nBitErr_erased + round(payloadLenBits * 0.5);
                    nPktErr  = nPktErr  + 1;
                    nErasure = nErasure + 1;

                elseif ~crc_ok
                    % --- CRC FAIL: данные есть, но ошибки присутствуют ---
                    % BER: сравниваем реальные биты
                    nBitErr_real = nBitErr_real + sum(bits_tx ~= bits_rx);
                    nPktErr  = nPktErr  + 1;
                    nCrcFail = nCrcFail + 1;

                else
                    % --- SUCCESS: пакет принят ---
                    % BER: сравниваем для полноты (обычно 0 ошибок)
                    nBitErr_real = nBitErr_real + sum(bits_tx ~= bits_rx);
                    nSuccess     = nSuccess + 1;
                end
            end

            % Суммарный BER по всем пакетам
            totalBits = Npkts * payloadLenBits;
            ber = (nBitErr_real + nBitErr_erased) / totalBits;
            per = nPktErr / Npkts;

            % Детальная статистика
            stats.nErasure    = nErasure;
            stats.nCrcFail    = nCrcFail;
            stats.nSuccess    = nSuccess;
            stats.nPktErr     = nPktErr;
            stats.perErasure  = nErasure  / Npkts;
            stats.perCrcFail  = nCrcFail  / Npkts;
            stats.perSuccess  = nSuccess  / Npkts;
            stats.ber         = ber;
            stats.per         = per;
            stats.Npkts       = Npkts;
            stats.payloadBits = payloadLenBits;
        end

        function printStats(~, stats)
            % printStats — вывод детальной статистики одного прогона.
            fprintf('  Пакетов всего : %d\n',    stats.Npkts);
            fprintf('  Успешных      : %d (%.1f%%)\n', stats.nSuccess,  stats.perSuccess  * 100);
            fprintf('  CRC fail      : %d (%.1f%%)\n', stats.nCrcFail,  stats.perCrcFail  * 100);
            fprintf('  Erasure       : %d (%.1f%%)\n', stats.nErasure,  stats.perErasure  * 100);
            fprintf('  PER           : %.4f\n',  stats.per);
            fprintf('  BER           : %.3e\n',  stats.ber);
        end
    end
end