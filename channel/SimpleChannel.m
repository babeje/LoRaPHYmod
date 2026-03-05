classdef SimpleChannel < handle
    properties
        fs       % частота дискретизации сигнала
        snr_dB   % SNR в дБ
        cfo_Hz   % частотный сдвиг (CFO)
    end
    
    methods
        function obj = SimpleChannel(fs, snr_dB, cfo_Hz)
            obj.fs     = fs;
            obj.snr_dB = snr_dB;
            obj.cfo_Hz = cfo_Hz;
        end
        
        function y = pass(obj, x)
            % x: комплексный baseband сигнал (столбец)
            if size(x,2) > 1
                x = x(:);
            end
            
            % 1) CFO
            t = (0:length(x)-1).' / obj.fs;
            x_cfo = x .* exp(1j*2*pi*obj.cfo_Hz*t);
            
            % 2) АБГШ
            y = obj.addAwgn(x_cfo, obj.snr_dB);
        end
    end
    
    methods (Access = protected)
        function y = addAwgn(~, x, snr_dB)
            % Реализация AWGN без Communications Toolbox
            snr_lin = 10^(snr_dB/10);
            P_sig   = mean(abs(x).^2);
            P_noise = P_sig / snr_lin;
            noise   = sqrt(P_noise/2) * (randn(size(x)) + 1j*randn(size(x)));
            y = x + noise;
        end
    end
end
