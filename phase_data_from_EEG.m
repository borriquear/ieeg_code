function [ phase_synchronization_matrix, euler_phase_differences_vector ] = phase_data_from_EEG( EEG_study )
%phase_data_vector return a vector with the phase difference between all
%pairs for a EEG struture
% Jaime Gomez Ramirez, 2016
%   output  phase_data_vector
%   input EEG data structure
%% load the vector with channel names
%  EEG_cut_BL_HYP.chanlocs.labels
% vector with the name of the channels names of the channels you want to compute connectivity between
EEG = EEG_study;
channel_labels = zeros(EEG_study.nbchan);
channel_labels = {EEG_study.chanlocs.labels};

% create complex Morlet wavelet
center_freq_bands = [2.5 5 10 20 40] %delta, theta, alpha,beta, gamma 
freq_index = 5; %1..5
center_freq =center_freq_bands(freq_index); %5 
time      = -1:1/EEG.srate:1;
wavelet   = exp(2*1i*pi*center_freq.*time) .* exp(-time.^2./(2*(4/(2*pi*center_freq))^2));
half_wavN = (length(time)-1)/2;

% FFT parameters
n_wavelet = length(time);
n_data    = EEG.pnts;
n_conv    = n_wavelet+n_data-1;

% FFT of wavelet
waveletX = fft(wavelet,n_conv);
waveletX = waveletX ./ max(waveletX);

% initialize output time-frequency data
phase_data = zeros(2,EEG.pnts);
real_data  = zeros(2,EEG.pnts);
%calculate the phase data vector
initchan = 2;
disp([ 'Calculating the phase for All channels ....' ])
for i=initchan:EEG_study.nbchan
    chan2idx = find(strcmpi(channel_labels(i),{EEG.chanlocs.labels}));
    % analytic signal of channel 2
    fft_data = fft(squeeze(EEG.data(chan2idx,:,1)),n_conv);
    as = ifft(waveletX.*fft_data,n_conv);
    as = as(half_wavN+1:end-half_wavN);
    % collect real and phase data
    phase_data(i-1,:) = angle(as);
    real_data(i-1,:)  = real(as);
end 
disp([ 'OK' ])

% phase angle differences
% loop to build the phase synchronization matrix for all pair of electrodes
phase_synchronization_matrix = zeros(EEG_study.nbchan-1, EEG_study.nbchan-1);
disp([ 'Calculating Synchronization between All channels ....' ])
index_t = 1;
[one twos] = size(EEG_study.times)
euler_phase_differences_vector = [];
for i=initchan:EEG_study.nbchan
    %indexphasediff = i;
    for j=i:EEG_study.nbchan
        phase_angle_differences = phase_data(i-1,:)- phase_data(j-1,:);
        % euler representation of angles
        euler_phase_differences = exp(1i*phase_angle_differences);

        % mean vector (in complex space)
        mean_complex_vector = mean(euler_phase_differences);
        euler_phase_differences_vector(index_t) = mean_complex_vector;
        index_t = index_t +1;
        % length of mean vector (this is the "M" from Me^ik, and is the measure of phase synchronization)
        phase_synchronization_matrix(i-1,j-1) = abs(mean_complex_vector);
        phase_synchronization_matrix(j-1,i-1) = abs(mean_complex_vector);
        disp([ 'Synchronization between ' ,channel_labels(i) ' and ' channel_labels(j)  ' is ' num2str(phase_synchronization_matrix(i-1,j-1)) ])
    end
end
disp([ 'OK' ])
%
%subplot(221)
%ispc_time =  abs(mean(euler_phase_differences_vector))
 

patient = 'BS27';
patient = 'CJ28';
patient = 'WSA30';
patient = 'SM31';
w_freqency = center_freq; %wavelet frequency
condition = 'EEG\_BL\_EO\_PRE';

condition = 'EEG\_BL\_HYP';
condition = 'EEG\_BL\_EC\_PRE';
condition = 'EEG\_BL\_EC\_POST';

typeimplant = 'Grid64 + Hem-strips';
typeimplant = 'Bitemp + burr hole';
typeimplant = 'bitmp+frontal+interhem-';
figure
polar_abs =  abs(mean(euler_phase_differences_vector));
polar_angle = mean(mean(phase_synchronization_matrix))
%subplot(221)
h=polar([0 polar_angle],[0 polar_abs]);
set(h,'linewidth',6,'color','g')
xlabel(['condition:' condition  ',  patient:' patient])
t=({[' Mean polar abs=' num2str(polar_abs) ', Mean polar angle =' num2str(polar_angle) ' Freq ' num2str(center_freq)]});
title(t);

figure
disp([ 'Drawingcolor matrix ....' ])

%t=({['ISPC, patient' num2str(x)];['text_B=' num2str(y)];['text_C=' num2str(y)]});
clims = [0 1];
imagesc(phase_synchronization_matrix,clims);
colorbar;
t=({['ISPC, patient:' patient '  Wavelet frequency=' num2str(w_freqency)];['Condition:' condition ' Implant:' typeimplant]});
title(t);
xlabel('Channel') % x-axis label
ylabel('Channel') % y-axis label});
title(t);
xlabel({[ ' Mean polar abs=' num2str(polar_abs) ', Mean polar angle =' num2str(polar_angle)]}) % x-axis label
ylabel('Channel') % y-axis label
disp({[ ' Mean polar abs=' num2str(polar_abs) ', Mean polar angle =' num2str(polar_angle)]})
%draw_chart_polar (patient, condition)
end
function [] = draw_chart_polar (patient, condition)
%tableimplnt 
%patient   Implant     condition Frquency     rho      angle(mean matrix)     
% 27     bitmp-36       HYP      2.5       0.15294   0.2275
%                                 5        0.22021   0.21205
%                                 10       0.32556   0.3283
%                                 20       0.22194   0.2059
%                                 40       0.23031   0.20992
% 31  bitmp+frtl+interh-104  HYP    2.5    0.041419   0.11016
%                                     5    0.04984    0.08506
%                                     10   0.042007   0.059913
%                                     20   0.051835   0.056487
%                                     40   0.051213   0.049183
% 31  bitmp+frtl+interh-104  EC_PRE   2.5  0.040331   0.11043
%                                       5  0.04106    0.074815
%                                      10  0.042311   0.055846
%                                      20  0.073293    0.076385
%                                      40  0.085356    0.078951
% 31                         EC_POST   2.5 0.042337    0.10519
%                                       5  0.041645     0.075323
%                                       10  0.042887 0.059411
%                                       20   0.075102 0.079908
%                                       40 0.080931  0.075033

%28   Grid-interh-104           HYP    2.5 0.10293     0.14029
%                                        5 0.1184      0.18012
%                                        10 0.11736     0.16295
%                                        20 0.086211     0.12211
%                                        40   0.069191   0.076082
%                               EC_POST    2.5 0.10211 0.1356
%                                          5   0.12179 0.19355
%                                          10  0.11029 0.16406
%                                          20  0.085205 0.11581
%                                          40  0.075528  0.079613
%WSA 30 bitemp-36                  HYP     2.5  0.21081  0.20657
%                                           5   0.21207   0.1995
%                                           10  0.15316   0.13953
%                                            20 0.16363   0.14275
%                                            40 0.25352   0.23426
%                                   EC_PRE   2.5 0.12669  0.14758
%                                            5  0.13687  0.14105
%                                            10 0.11378    0.12133
%                                            20  0.16442    0.1461
%                                              40 0.26786 0.24863
%                                     EC_POST  2.5 0.21329 0.21079
%                                                 5 0.21943 0.209
%                                                 10 0.15226  0.13852
%                                                 20 0.1574  0.13708
%                                                   40 0.25884 0.23983

% 

patient = 'CJ28 (High-Hyp)'
patient = 'WSA30 (High-Hyp)'
patient = 'SM31 (Low-Hyp)'
sm31 = findstr(patient,  'SM31 (Low-Hyp)')
cj28 = findstr(patient, 'CJ28 (High-Hyp)')
wsa30 = findstr(patient, 'WSA30 (High-Hyp)')
if sm31 == 1
    typeimplant = 'bitmp+frontal+interhem-104e'
    rhos_ec_pre = [0.040331 0.04106  0.042311 0.073293 0.085356]
    angles_ec_pre = [0.11043 0.074815 0.055846 0.076385  0.078951]
    angles_ec_post = [0.10519 0.075323 0.059411 0.079908 0.075033]
    rhos_hyp = [0.041419 0.04984 0.042007 0.051835  0.051213]
    angles_hyp = [0.11016 0.08506  0.059913 0.056487 0.049183]
    rhos = [rhos_ec_pre;rhos_hyp]
    angles = [angles_ec_pre(1) angles_ec_post(1) angles_hyp(1);angles_ec_pre(2) angles_ec_post(2) angles_hyp(2);angles_ec_pre(3) angles_ec_post(3) angles_hyp(3);angles_ec_pre(4) angles_ec_post(4) angles_hyp(4);angles_ec_pre(5) angles_ec_post(5) angles_hyp(5)]
    bar(angles)
    ylabel('Phase synchronization')
    set(gca,'xticklabel',{'delta','theta','alpha','beta','gamma'})
    legend('EC\_PRE','EC\_POST','HYP')
    title({['Frequency-based syOnchronization. Patient:' patient ] ['Implant:' typeimplant]})
end
if cj28 == 1 
    typeimplant = 'Grid+interhem-88e'
    rhos_ec_post = [];
    angles_ec_post = [ 0.1356 0.19355 0.16406 0.11581 0.079613];
    rhos_hyp = [];
    angles_hyp = [0.14029  0.18012 0.16295  0.12211 0.076082];
    angles = [angles_ec_post(1) angles_hyp(1);angles_ec_post(2) angles_hyp(2);angles_ec_post(3) angles_hyp(3);angles_ec_post(4) angles_hyp(4);angles_ec_post(5) angles_hyp(5)]
    bar(angles);
    ylabel('Phase synchronization');
    set(gca,'xticklabel',{'delta','theta','alpha','beta','gamma'});
    legend('EC\_POST','HYP');
    title({['Frequency-based synchronization. Patient:' patient ] ['Implant:' typeimplant]});
end
if wsa30 ==1
    typeimplant = 'bitmp-36e'
    rhos_ec_post = [];
    rhos_hyp = [];
    angles_hyp = [  0.20657  0.1995 0.13953  0.14275 0.23426 ];
    angles_ec_post = [ 0.21079  0.209  0.13852 0.13708 0.23983];
    angles_ec_pre = [ 0.14758 0.14105 0.12133  0.1461 0.24863];
    angles = [angles_ec_pre(1) angles_ec_post(1) angles_hyp(1);angles_ec_pre(2) angles_ec_post(2) angles_hyp(2);angles_ec_pre(3) angles_ec_post(3) angles_hyp(3);angles_ec_pre(4) angles_ec_post(4) angles_hyp(4);angles_ec_pre(5) angles_ec_post(5) angles_hyp(5)]
    angles_ec_pre= []
    bar(angles);
    ylabel('Phase synchronization');
    set(gca,'xticklabel',{'delta','theta','alpha','beta','gamma'});
    legend('EC\_PRE', 'EC\_POST','HYP');
    title({['Frequency-based synchronization. Patient:' patient ] ['Implant:' typeimplant]});
end
end
