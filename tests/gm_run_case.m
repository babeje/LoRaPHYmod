function res = gm_run_case(cs)
% gm_run_case  Детерминированное исполнение одного Golden-Master-кейса.
%
% Назначение:
%   Исполняет вычислительное ядро проекта (LoRaModem + канал + LoRaSimulator
%   либо LoRaModem.estimateChannel) с фиксированным seed и возвращает
%   числовые результаты, пригодные для строгой регрессионной сверки.
%
%   Обвязка намеренно гарантирует КОМПОНЕНТЫ ядра, а не тела сценариев
%   (циклы/графики/сохранение), поскольку именно компоненты затрагиваются
%   рефакторингом. Логика сценариев (run_snr_sweep, run_channel_estimation_demo)
%   подлежит выносу и не является предметом фиксации Golden Master.
%
% Детерминизм:
%   rng(cs.seed, 'twister') задаёт единый глобальный поток на весь кейс.
%   Все стохастические потребители (генерация бит в LoRaSimulator.run,
%   реализации замираний RayleighTDLChannel при Seed=[], аддитивный шум)
%   читают этот поток последовательно ⇒ при неизменном порядке вызовов
%   результат воспроизводим бит-в-бит.
%
% Вход:
%   cs [struct] — спецификация кейса (см. gm_define_cases в run_golden_master.m).
%                 Поле cs.kind ∈ {'link', 'chanest'} выбирает ветвь.
%
% Выход:
%   res [struct] — числовые массивы результатов (состав зависит от cs.kind).
%
% Зависимости: core/LoRaModem.m, core/LoRaSimulator.m,
%              channel/RayleighTDLChannel.m, channel/AwgnChannel.m
% Совместимость: MATLAB R2024a, без дополнительных тулбоксов
%                (при cs.FastMode = true).
% ---------------------------------------------------------------

    rng(cs.seed, 'twister');   % Генерация детерменированного 
                               % случайного потока чисел
                               % для воспроизводимости

    switch cs.kind
        case 'link'
            res = local_run_link(cs);
        case 'chanest'
            res = local_run_chanest(cs);
        otherwise
            error('gm_run_case:unknownKind', ...
                'Неизвестный тип кейса: %s', cs.kind);
    end
end

% ------------------------------------------------------------------
%  Ветвь 'link' — свип BER/PER по (скорость × ОСШ) через полный тракт
% ------------------------------------------------------------------
function res = local_run_link(cs)
    c  = 299792458;                 % скорость света, м/с
    nV = numel(cs.v_mps_list);
    nS = numel(cs.snr_list);

    BER      = zeros(nV, nS);
    PER      = zeros(nV, nS);
    nErasure = zeros(nV, nS);       % сырые счётчики исходов — целочисленные,
    nCrcFail = zeros(nV, nS);       % воспроизводимы точно (основа строгой сверки)
    nSuccess = zeros(nV, nS);
    nPktErr  = zeros(nV, nS);

    for vi = 1:nV
        fd_Hz = (cs.v_mps_list(vi) / c) * cs.rf_freq;   % допплеровское расширение

        for si = 1:nS
            snr_dB = cs.snr_list(si);

            modem = LoRaModem(cs.rf_freq, cs.sf, cs.bw, cs.fs, ...
                'CR',          cs.CR, ...
                'HasHeader',   true, ...
                'UseCRC',      true, ...
                'PreambleLen', cs.PreambleLen, ...
                'FastMode',    cs.FastMode);

            switch cs.channel
                case 'tdl'
                    ch = RayleighTDLChannel(cs.fs, snr_dB, 0, ...
                        'SF',         cs.sf, ...
                        'BW',         cs.bw, ...
                        'DopplerHz',  fd_Hz, ...
                        'PathDelays', cs.tdlDelays, ...
                        'PathGains',  cs.tdlGains, ...
                        'Seed',       []);           % реализации из глобального потока
                case 'awgn'
                    ch = AwgnChannel(cs.fs, snr_dB);
                otherwise
                    error('gm_run_case:badChannel', ...
                        'Неизвестный канал: %s', cs.channel);
            end

            sim = LoRaSimulator(modem, ch);
            [ber, per, st] = sim.run(cs.Npkts, cs.payloadBits);

            BER(vi, si)      = ber;
            PER(vi, si)      = per;
            nErasure(vi, si) = st.nErasure;
            nCrcFail(vi, si) = st.nCrcFail;
            nSuccess(vi, si) = st.nSuccess;
            nPktErr(vi, si)  = st.nPktErr;
        end
    end

    res = struct('BER', BER, 'PER', PER, ...
                 'nErasure', nErasure, 'nCrcFail', nCrcFail, ...
                 'nSuccess', nSuccess, 'nPktErr', nPktErr);
end

% ------------------------------------------------------------------
%  Ветвь 'chanest' — Монте-Карло качества ML/MMSE-оценки канала
% ------------------------------------------------------------------
function res = local_run_chanest(cs)
    os = cs.fs / cs.bw;
    Ns = 2^cs.sf * os;              % отсчётов на символ

    % Эталонный TX-сигнал (детерминирован; payload не влияет на оценку h)
    modem_gen = LoRaModem(cs.rf_freq, cs.sf, cs.bw, cs.fs, ...
        'CR', cs.CR, 'HasHeader', true, 'UseCRC', true, ...
        'PreambleLen', 8, 'FastMode', true);

    bits_tx        = logical(randi([0 1], cs.payloadBits, 1));
    [txSig, ~, ~]  = modem_gen.modulate(bits_tx);
    s_ref          = txSig(1:Ns);
    Ntx            = numel(txSig);

    nS = numel(cs.snr_list);
    nC = numel(cs.preamble_configs);

    rmse_amp = zeros(nS, nC);
    var_cplx = zeros(nS, nC);

    for si = 1:nS
        snr_dB  = cs.snr_list(si);
        snr_lin = 10^(snr_dB / 10);
        nvar    = 1 / snr_lin;      % дисперсия шума на комплексный отсчёт

        for ci = 1:nC
            Np = cs.preamble_configs(ci);
            modem_est = LoRaModem(cs.rf_freq, cs.sf, cs.bw, cs.fs, ...
                'CR', cs.CR, 'HasHeader', true, 'UseCRC', true, ...
                'PreambleLen', Np, 'FastMode', true);

            h_est_arr  = complex(zeros(cs.N_mc, 1));
            h_true_arr = complex(zeros(cs.N_mc, 1));

            for t = 1:cs.N_mc
                h_true = (randn + 1j*randn) / sqrt(2);           % CN(0,1)
                noise  = sqrt(nvar/2) * (randn(Ntx,1) + 1j*randn(Ntx,1));
                rxSig  = h_true * txSig + noise;                 % плоский канал
                [h_est, ~] = modem_est.estimateChannel(rxSig, snr_dB, s_ref);
                h_est_arr(t)  = h_est;
                h_true_arr(t) = h_true;
            end

            rmse_amp(si, ci) = sqrt(mean((abs(h_est_arr) - abs(h_true_arr)).^2));
            var_cplx(si, ci) = mean(abs(h_est_arr - h_true_arr).^2);
        end
    end

    res = struct('rmse_amp', rmse_amp, 'var_cplx', var_cplx);
end
