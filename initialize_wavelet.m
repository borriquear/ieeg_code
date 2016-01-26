function [half_of_wavelet_size, n_wavelet, n_data, n_convolution, wavelet_cycles] =   initialize_wavelet(time, EEG_study)
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters
n_wavelet     = length(time);
n_data        = EEG_study.pnts*EEG_study.trials;
n_convolution = n_wavelet+n_data-1;
wavelet_cycles= 4.5;
end