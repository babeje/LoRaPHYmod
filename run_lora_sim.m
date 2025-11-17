clear; clc;

% Параметры PHY
rf_freq = 868e6;
sf      = 7;
bw      = 125e3;
fs      = 1e6;

% Параметры канала
snr_dB  = -10;   % SNR в дБ
cfo_Hz  = 100;   % CFO, Гц

% Длина полезной нагрузки (в битах)
payloadLenBits = 128;

% Кол-во пакетов
Npkts = 200;

% Создаём модем и канал
modem   = LoRaModem(rf_freq, sf, bw, fs, ...
                    'CR', 4, ...
                    'HasHeader', true, ...
                    'UseCRC', true, ...
                    'PreambleLen', 8, ...
                    'FastMode', false);

channel = SimpleChannel(fs, snr_dB, cfo_Hz);

% Имитационная модель
sim     = LoRaSimulator(modem, channel);

[ber, per] = sim.run(Npkts, payloadLenBits);

fprintf('SNR = %.1f dB, CFO = %.1f Hz\n', snr_dB, cfo_Hz);
fprintf('BER = %.3e, PER = %.3f\n', ber, per);
