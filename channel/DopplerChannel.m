classdef DopplerChannel < SimpleChannel
    % DopplerChannel: AWGN + (CFO + Doppler)
    % Doppler: f_d = (v/c)*fc*cos(theta)
    % Можно задать либо напрямую v/theta, либо напрямую fd_Hz.

    properties
        fc_Hz = 868e6            % несущая, Гц
        v_mps = 0                % скорость, м/с
        theta_rad = 0            % угол
        doppler_rate_Hzps = 0    % изменение допплера во времени, Гц/с (0 = постоянный)
    end

    methods
        function obj = DopplerChannel(fs, snr_dB, cfo_Hz, fc_Hz, v_mps, theta_rad, doppler_rate_Hzps)
            obj@SimpleChannel(fs, snr_dB, cfo_Hz);

            if nargin >= 4 && ~isempty(fc_Hz), obj.fc_Hz = fc_Hz; end
            if nargin >= 5 && ~isempty(v_mps), obj.v_mps = v_mps; end
            if nargin >= 6 && ~isempty(theta_rad), obj.theta_rad = theta_rad; end
            if nargin >= 7 && ~isempty(doppler_rate_Hzps), obj.doppler_rate_Hzps = doppler_rate_Hzps; end
        end

        function y = pass(obj, x)
            if size(x,2) > 1
                x = x(:);
            end

            t = (0:length(x)-1).' / obj.fs;

            % Допплер (постоянная составляющая)
            c = 299792458; % м/с
            fd0 = (obj.v_mps / c) * obj.fc_Hz * cos(obj.theta_rad);

            % Дополнительный доплеровский эффект с изменяющимися во времени параметрами:
            % fd(t) = fd0 + rate*t
            fd_t = fd0 + obj.doppler_rate_Hzps .* t;

            % Общий сдвиг = CFO + Doppler(t)
            f_total = obj.cfo_Hz + fd_t;

            % Фаза = 2*pi*∫f_total dt ≈ 2*pi*cumsum(f_total)/fs
            phase = 2*pi*cumsum(f_total) / obj.fs;

            x_shifted = x .* exp(1j * phase);

            % АБГШ
            y = obj.addAwgn(x_shifted, obj.snr_dB);
        end
    end
end
