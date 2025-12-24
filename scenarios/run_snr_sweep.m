clear; clc;

projectRoot = fileparts(mfilename('fullpath'));
projectRoot = fileparts(projectRoot);
addpath(genpath(projectRoot));

% параметры PHY
rf_freq = 868e6;
sf      = 10;
bw      = 125e3;
fs      = 1e6;

% Параметры модема
CR          = 4;
HasHeader   = true;
UseCRC      = true;
PreambleLen = 8;
FastMode    = true;   % можно true для ускорения

% Симуляция
Npkts          = 200;
payloadLenBits = 256;

% Настройки разброса по ОСШ
snr_list = -10:2:2;

% Настройки доплера
cases(1).name = "v=0";
cases(1).v_mps = 0;
cases(1).theta = 0;
cases(1).fd_rate = 0;

cases(2).name = "v=30, fd_rate=5000";
cases(2).v_mps = 30;
cases(2).theta = 0;
cases(2).fd_rate = 2500;

% Можно добавить еще сценариев
% cases(3).name = "v=30, fd_rate=20000";
% cases(3).v_mps = 30; cases(3).theta = 0; cases(3).fd_rate = 20000;

% Результаты
PER = zeros(numel(cases), numel(snr_list));
BER = zeros(numel(cases), numel(snr_list));

for ci = 1:numel(cases)
    fprintf("\\n=== Case: %s ===\\n", cases(ci).name);

    for si = 1:numel(snr_list)
        snr_dB = snr_list(si);
        cfo_Hz = 0;

        modem = LoRaModem(rf_freq, sf, bw, fs, ...
            'CR', CR, 'HasHeader', HasHeader, 'UseCRC', UseCRC, ...
            'PreambleLen', PreambleLen, 'FastMode', FastMode);

        channel = DopplerChannel(fs, snr_dB, cfo_Hz, rf_freq, ...
            cases(ci).v_mps, cases(ci).theta, cases(ci).fd_rate);

        sim = LoRaSimulator(modem, channel);
        [ber, per] = sim.run(Npkts, payloadLenBits);

        BER(ci, si) = ber;
        PER(ci, si) = per;

        fprintf("SNR=%5.1f dB  BER=%9.3e  PER=%6.3f\\n", snr_dB, ber, per);
    end
end

figure;
plot(snr_list, PER(1,:), '-o'); hold on;
for ci = 2:numel(cases)
    plot(snr_list, PER(ci,:), '-o');
end
grid on;
xlabel('SNR (dB)');
ylabel('PER');
legend(string({cases.name}), 'Location', 'best');
title(sprintf('PER, SF=%d, BW=%.0f kHz, Количество пакетов=%d, Длина полезной нагрузки=%d bits', sf, bw/1e3, Npkts, payloadLenBits));

outDir = fullfile(projectRoot, 'results', 'data');
if ~exist(outDir), mkdir(outDir); end

ts = datestr(now, 'yyyymmdd_HHMMSS');
save(fullfile(outDir, "snr_sweep_" + ts + ".mat"), 'snr_list', 'cases', 'BER', 'PER', ...
     'rf_freq','sf','bw','fs','Npkts','payloadLenBits','CR','HasHeader','UseCRC','PreambleLen','FastMode');
