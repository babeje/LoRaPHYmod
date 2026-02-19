classdef MultipathChannel < SimpleChannel
    % MultipathChannel: частотно-селективный релеевский канал с медленными замираниями
    % Добавляет многолучевость, доплеровское расширение, а затем CFO и AWGN.
    
    properties
        rayleighChan        % объект comm.RayleighChannel
        pathDelays          % задержки лучей (сек)
        pathGains           % относительные мощности лучей (дБ)
        maxDopplerShift     % максимальный доплеровский сдвиг (Гц)
    end
    
    methods
        function obj = MultipathChannel(fs, snr_dB, cfo_Hz, varargin)
            % Конструктор
            % Обязательные: fs, snr_dB, cfo_Hz
            % Дополнительные параметры (в виде пар 'Name',Value):
            %   'PathDelays'       - вектор задержек (сек), по умолчанию [0, 15e-6]
            %   'PathGains'        - вектор мощностей (дБ), по умолчанию [0, 0]
            %   'MaxDopplerShift'  - максимальный доплер (Гц), по умолчанию 87 (для v=30 м/с)
            %   'Seed'             - начальное число для генератора (по умолчанию случайное)
            
            obj@SimpleChannel(fs, snr_dB, cfo_Hz);
            
            % Параметры по умолчанию (частотно-селективный канал с Bc ~ 56.5 кГц)
            defaultPathDelays = [0, 15e-6];      % задержки 0 и 15 мкс
            defaultPathGains = [0, 0];            % равные мощности
            defaultMaxDoppler = 87;                % fm = v * fc / c, v=30 м/с, fc=868e6
            defaultSeed = randi(2^32-1);
            
            p = inputParser;
            addParameter(p, 'PathDelays', defaultPathDelays);
            addParameter(p, 'PathGains', defaultPathGains);
            addParameter(p, 'MaxDopplerShift', defaultMaxDoppler);
            addParameter(p, 'Seed', defaultSeed);
            parse(p, varargin{:});
            
            obj.pathDelays = p.Results.PathDelays;
            obj.pathGains = p.Results.PathGains;
            obj.maxDopplerShift = p.Results.MaxDopplerShift;
            
            % Создаем объект многолучевого канала
            obj.rayleighChan = comm.RayleighChannel(...
                'SampleRate', obj.fs, ...
                'PathDelays', obj.pathDelays, ...
                'AveragePathGains', obj.pathGains, ...
                'MaximumDopplerShift', obj.maxDopplerShift, ...
                'DopplerSpectrum', doppler('Jakes'), ...
                'NormalizePathGains', true, ...
                'RandomStream', 'mt19937ar with seed', ...
                'Seed', p.Results.Seed);
        end
        
        function y = pass(obj, x)
            % Переопределение метода pass: многолучевость → CFO → AWGN
            if size(x,2) > 1
                x = x(:);
            end
            
            % Многолучевое распространение с замираниями
            x_mp = obj.rayleighChan(x);
            
            % Добавляем CFO (как в SimpleChannel)
            t = (0:length(x_mp)-1).' / obj.fs;
            x_cfo = x_mp .* exp(1j*2*pi*obj.cfo_Hz*t);
            
            % Добавляем AWGN
            y = obj.addAwgn(x_cfo, obj.snr_dB);
        end
    end
end