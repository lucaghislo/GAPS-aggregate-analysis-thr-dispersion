function [Thr_eff_keV] = conversion_DACthru_keV(G_ch_ADCu_DACinju, Thr_DAC, type)
% Steps:
% Thr_eff_mV = Thr_max_mV - Thr_DAC * conv_mV_DACthru [mv]
% G_ch_mV_keV = G_ch_ADCu_DACinju * conv_mV_ADCu / conv_keV_DACinju [mV/keV]
% G_sh unit is 1 or ADCu/mV?
% Thr_eff_keV = Thr_eff_mV * G_sh / G_ch_mV_keV [keV]

% conv_keV_DACinju = 2.35*.8297; % [kev / DAC_inj code]; (2.35 is for width at half-height)
% conv_keV_DACinju = 2.35*.838; % [keV / DAC_inj code] from ENC.m, our gamma is = 0.885
conv_keV_DACinju = .841; % [keV / DAC_inj code]

FF_conv_mV_DACthru = 1.36; % [mV/DAC_thru]
TT_conv_mV_DACthru = 1.58; % [mV/DAC_thru] % corretto da 1.20 a 1.58
SS_conv_mV_DACthru = 1.03; % [mV/DAC_thru]
FF_max_mV = 320; % [mV]
TT_max_mV = 280; % [mV]
SS_max_mV = 240; % [mV]

if strcmp(type, 'FF')
    max_mV = FF_max_mV;
    conv_mV_DACthru = FF_conv_mV_DACthru;
elseif strcmp(type, 'TT')
    max_mV = TT_max_mV;
    conv_mV_DACthru = TT_conv_mV_DACthru;
else
    max_mV = SS_max_mV;
    conv_mV_DACthru = SS_conv_mV_DACthru;
end

G_sh = 5.14; % [ADCu/mV]
conv_mV_ADCu = 1.72; % [mV/ADCu]

Thr_eff_mV = max_mV - Thr_DAC * conv_mV_DACthru;
G_ch_mV_keV = G_ch_ADCu_DACinju * conv_mV_ADCu / conv_keV_DACinju;
Thr_eff_keV = Thr_eff_mV * G_sh / G_ch_mV_keV;
end

