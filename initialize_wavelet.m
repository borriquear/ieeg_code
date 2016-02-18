function [half_of_wavelet_size, n_wavelet, n_data, n_convolution, wavelet_cycles, time] =   initialize_wavelet(EEG_study)
% initialize_wavelet  initialize complex wavelet parameters
% half_of_wavelet_size = (length(time)-1)/2;
% % FFT parameters
% time = -1:1/EEG_study.srate:1; time duration of the complext wavelet, it has
%   to have same sampling rate that the EEG signal
% n_wavelet     = length(time);
% n_data        = EEG_study.pnts*EEG_study.trials;
% n_convolution = n_wavelet+n_data-1;
% wavelet_cycles= 4.5;
time = -1:1/EEG_study.srate:1;
half_of_wavelet_size = (length(time)-1)/2;
% FFT parameters
n_wavelet     = length(time);
n_data        = EEG_study.pnts*EEG_study.trials;
n_convolution = n_wavelet+n_data-1;
wavelet_cycles= 4.5;
end