classdef AwgnChannel < SimpleChannel
    methods
        function obj = AwgnChannel(fs, snr_dB)
            obj@SimpleChannel(fs, snr_dB, 0); % CFO=0
        end
    end
end