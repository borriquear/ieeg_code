%% phase-based connectivity analysis
% Calculate ISPC over time, for one single tril in long epochs, e.g.,
% reting state
%load the data set
%load 'C:\Users\Jaime\Documents\BIAL PROJECT\patients\Cordeiro, Jessica\Session 1 day 1 20 October\EEG_cut_BL_EC_POST_cj28_10202015_s1.mat'
%load 'C:\Users\Jaime\Documents\BIAL PROJECT\patients\Button, Steven\Session 2 22 October\EEG_cut_BL_HYP_bs27_10222015_s2.mat'
%load 'C:\Users\Jaime\Documents\BIAL PROJECT\patients\twh31\EEG_cut_BL_EC_PRE_sm31_12012015_s1.mat'
load 'C:\Users\Jaime\Documents\BIAL PROJECT\patients\twh31\EEG_cut_BL_EC_POST_sm31_12012015_s1.mat'

%load sampleEEGdata
%EEG_study = EEG_cut_BL_HYP;
%EEG_study = EEG_cut_BL_EC_POST;
%EEG_study = EEG_cut_BL_EC_PRE;
EEG_study = EEG_cut_BL_EC_POST;
%%
% phase angle differences
% loop to build the phase synchronization matrix for all pair of electrodes
[phase_synchronization_matrix, euler_phase_differences_vector ]  = phase_data_from_EEG(EEG_study);
return;
% phase_angle_differences = phase_data(2,:)-phase_data(1,:);
% 
% % euler representation of angles
% euler_phase_differences = exp(1i*phase_angle_differences);
% 
% % mean vector (in complex space)
% mean_complex_vector = mean(euler_phase_differences);
% 
% % length of mean vector (this is the "M" from Me^ik, and is the measure of phase synchronization)
% phase_synchronization = abs(mean_complex_vector);
% 
% disp([ 'Synchronization between ' channel1 ' and ' channel2 ' is ' num2str(phase_synchronization) '!' ])
% 
% % of course, this could all be done on one line:
% %phase_synchronization = abs(mean(exp(1i*(phase_data(2,:)-phase_data(1,:)))));
% %update [hasesynchronization matrix
% phase_synchronization_matix(i,j) = phase_synchronization

% notice that the order of subtraction is meaningless (see below), which means that this measure of synchronization is non-directional!
phase_synchronization_backwards = abs(mean(exp(1i*(phase_data(1,:)-phase_data(2,:)))));


% % now plot mean vector
% subplot(221)
% hold on
% h=polar([0 angle(mean_complex_vector)],[0 phase_synchronization]);
% set(h,'linewidth',6,'color','g')
% 
% subplot(222)
% hold on
% h=polar([0 mean(mean(mean_complex_matrix_angle))],[0 mean(mean(phase_synchronization_matrix))]);
% set(h,'linewidth',6,'color','g')
% 
% %% phase clustering is phase-invariant
% 
% figure(2), clf
% subplot(221)
% polar(repmat(phase_data(2,:)-phase_data(1,:),1,2)',repmat([0 1],1,EEG.pnts)','k');
% title([ 'Phase synchronization: ' num2str(abs(mean(exp(1i*(diff(phase_data,1)))))) ])
% 
% new_phase_data = phase_data;
% for i=2:4
%     subplot(2,2,i)
%     
%     % add random phase offset
%     new_phase_data(1,:) = new_phase_data(1,:)+rand*pi;
%     
%     % plot again
%     polar(repmat(new_phase_data(2,:)-new_phase_data(1,:)+pi/2,1,2)',repmat([0 1],1,EEG.pnts)','k');
%     title([ 'Phase synchronization: ' num2str(abs(mean(exp(1i*(diff(new_phase_data,1)))))) ])
% end
% 
% %% 
% 
% %% ISPC over time vs. over trials
% 
% % FFT parameters
% n_wavelet = length(time);
% n_data    = EEG.pnts*EEG.trials;
% n_conv    = n_wavelet+n_data-1;
% 
% 
% % initialize output time-frequency data
% phase_data = zeros(2,EEG.pnts,EEG.trials);
% 
% 
% % FFT of wavelet (need to redo FFT because different n_conv)
% waveletX = fft(wavelet,n_conv);
% waveletX = waveletX ./ max(waveletX);
% 
% 
% % analytic signal of channel 1
% fft_data = fft(reshape(EEG.data(chan1idx,:,:),1,[]),n_conv);
% as = ifft(waveletX.*fft_data,n_conv);
% as = as(half_wavN+1:end-half_wavN);
% as = reshape(as,EEG.pnts,EEG.trials);
% 
% % collect real and phase data
% phase_data(1,:,:) = angle(as);
% 
% % analytic signal of channel 1
% fft_data = fft(reshape(EEG.data(chan2idx,:,:),1,[]),n_conv);
% as = ifft(waveletX.*fft_data,n_conv);
% as = as(half_wavN+1:end-half_wavN);
% as = reshape(as,EEG.pnts,EEG.trials);
% 
% % collect real and phase data
% phase_data(2,:,:) = angle(as);
% 
% 
% figure(3), clf
% 
% % ISPC over trials
% subplot(211)
% ispc_trials = abs(mean(exp(1i*diff(phase_data,[],1)
% plot(EEG.times,ispc_trials)
% set(gca,'xlim',[-200 1200])
% xlabel('Time (ms)'), ylabel('ISPC')
% 
% % ISPC over time
% subplot(212)
% ispc_time = squeeze( abs(mean(exp(1i*diff(phase_data,[],1)))))
% plot(1:EEG.trials,ispc_time)
% xlabel('Trials'), ylabel('ISPC')

%% end.
