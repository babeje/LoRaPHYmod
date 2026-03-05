clear; clc;

projectRoot = fileparts(mfilename('fullpath'));
projectRoot = fileparts(projectRoot);
addpath(genpath(projectRoot));

% Параметры PHY
rf_freq = 868e6;
sf      = 9;
bw      = 125e3;
fs      = 1e6;

% Параметры канала
snr_dB  = -8;   % SNR в дБ
cfo_Hz  = 0;   % CFO, Гц

% Длина полезной нагрузки (в битах)
payloadLenBits = 128;

% Кол-во пакетов
Npkts = 50;

% Модем и канал
modem   = LoRaModem(rf_freq, sf, bw, fs, ...
                    'CR', 4, ...
                    'HasHeader', true, ...
                    'UseCRC', true, ...
                    'PreambleLen', 8, ...
                    'FastMode', false);

v_mps     = 30;     % скорость, м/с
theta_rad = 0;      % 0 = летим на приёмник (макс допплер)
fd_rate   = 5000;      % пока без изменения допплера во времени

channel = DopplerChannel(fs, snr_dB, cfo_Hz, rf_freq, v_mps, theta_rad, fd_rate);

% Имитационная модель
sim     = LoRaSimulator(modem, channel);

[ber, per] = sim.run(Npkts, payloadLenBits);

fprintf('SNR = %.1f dB, CFO = %.1f Hz\n', snr_dB, cfo_Hz);
fprintf('BER = %.3e, PER = %.3f\n', ber, per);
