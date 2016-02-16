function [ output_args ] = phasebasedconnectivity( input_args )
%% phasebasedconnectivity  Phase based connectivity analysis
% Returns phase_synchronization_matrix, euler_phase_differences_vector
% 1. Load data
%    Call to initialize_EEG_variables() and initialize_wavelet(EEG_study)
% 2. Calculate phase for all channels, obtain phase_synchronization_matrix
    % Call phasedatafromEEG(EEG, half_of_wavelet_size, n_wavelet, n_data, n_convolution, wavelet_cycles, time );
% 
%% 1. Load epoch and Initialize data
disp('Loading data....')
[myfullname, EEG_study, channel_labels] = initialize_EEG_variables();
EEG = EEG_study;
[half_of_wavelet_size, n_wavelet, n_data, n_convolution, wavelet_cycles, time] = initialize_wavelet(EEG);
% draw_all_single_correlations = 0 do not draw all correlation chart,
% justmatrix, set to 0 when calculating the entire matrix
draw_all_single_correlations = 0;
times2plot =  dsearchn(EEG.times',[ EEG.times(1) EEG.times(end)]');
trial2plot = EEG.trials;
%centerfreq = 10;
[pathstr, subjandcond, matext] = fileparts(myfullname);
%To plot only one channel: initchan = chantouse
%To plot from channel 12 to 23, initchan = 12, chantouse = 23;
center_freq_bands = [2.5 5 10 20 40] %delta, theta, alpha,beta, gamma
initchan =  2;
chantouse = EEG.nbchan;
disp('Data Loaded OK!')

disp('Calculating the phase synchronization matrix ...')
[phase_synchronization_matrix, euler_phase_differences_vector ]  = phasedatafromEEG(EEG, center_freq_bands, half_of_wavelet_size, n_wavelet, n_data, n_convolution, wavelet_cycles, time,myfullname);
disp('DONE!');
end

function [phase_synchronization_matrix, euler_phase_differences_vector ]  = phasedatafromEEG(EEG, center_freq_bands, half_wavN, n_wavelet, n_data, n_conv, wavelet_cycles, time,myfullname)
%% phasedatafromEEG calculates the phase_synchronization_matrix and euler_phase_differences_vector
% for a EEG object
channel_labels = {EEG.chanlocs.labels};
% Frequencies to study

disp(['Calculating phase based synchrony for ' num2str(size(center_freq_bands,2)) ' frequency bands '] )
%
for freq_index=1:size(center_freq_bands,2)
    
    center_freq = center_freq_bands(freq_index); %5
    wavelet   = exp(2*1i*pi*center_freq.*time) .* exp(-time.^2./(2*(4/(2*pi*center_freq))^2));
    % FFT of wavelet
    waveletX = fft(wavelet,n_conv);
    waveletX = waveletX ./ max(waveletX);
    % initialize output time-frequency data
    phase_data = zeros(2,EEG.pnts);
    real_data  = zeros(2,EEG.pnts);
    %calculate the phase data vector
    initchan = 2;
    disp([ 'Calculating the phase for frequency = ' num2str(center_freq) ' for All channels ....' ])
    for i=initchan:EEG.nbchan
        chan2idx = find(strcmpi(channel_labels(i),{EEG.chanlocs.labels}));
        % analytic signal of channel 2
        fft_data = fft(squeeze(EEG.data(chan2idx,:,1)),n_conv);
        as = ifft(waveletX.*fft_data,n_conv);
        as = as(half_wavN+1:end-half_wavN);
        % collect real and phase data
        phase_data(i-1,:) = angle(as);
        real_data(i-1,:)  = real(as);
    end
    
    % phase angle differences
    % loop to build the phase synchronization matrix for all pair of electrodes
    phase_synchronization_matrix = zeros(EEG.nbchan-1, EEG.nbchan-1);
    disp([ 'Calculating Synchronization between All channels ....' ])
    index_t = 1;
    [one twos] = size(EEG.times)
    euler_phase_differences_vector = [];
    for i=initchan:EEG.nbchan
        %indexphasediff = i;
        for j=i:EEG.nbchan
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
    disp([ 'Synchronization between All channels  DONE!' ])
    % Get patient, condition from myfullname
    [pathstr, subjandcond, matext] = fileparts(myfullname);
    [BL remain] = strtok(subjandcond, 'EEG_cut_');
    [condition remain] = strtok(remain, '_');
    [patient remain]  = strtok(remain, '_');

     polar_abs =  abs(mean(euler_phase_differences_vector));
     polar_angle = mean(mean(phase_synchronization_matrix))
     % draw_polar_coords =1 if want to display the polar coordinates
     draw_polar_coords = 0
     if draw_polar_coords == 1
         figure
         h=polar([0 polar_angle],[0 polar_abs]);
         set(h,'linewidth',6,'color','g')
         xlabel(['condition:' condition  ',  patient:' patient])
         t=({[' Mean polar abs=' num2str(polar_abs) ', Mean polar angle =' num2str(polar_angle) ' Freq ' num2str(center_freq)]});
        title(t);
     end
     disp([' Saving phase synchronization matrix for Freq=' num2str(center_freq) ' in ' myfullname])
     %Save object  corrmat as mat file
     [mat_pathstr,mat_name,mat_ext] = fileparts(myfullname);
%     %typeofcorr =strcat(typeofcorr, '_');
%     ftypeofcorr =strcat(typeofcorr, '_Fq_');
%     %freq = sprintf('%d_', centerfreq);
%     ftypeofcorr =strcat(ftypeofcorr,num2str(centerfreq));
%     ftypeofcorr = strcat(ftypeofcorr,'_');
%     corr_mat_name = strcat(ftypeofcorr,mat_name);
%     %create directory if does not exist 
    cd(mat_pathstr);
    if ~exist(mat_name, 'dir')
        mkdir(mat_name);
    end
    mat_pathstr =  fullfile(mat_pathstr,mat_name);
    %mat_pathstr =  strcat(mat_pathstr,mat_name);
    cd(mat_pathstr);
    phaseandfreq = strcat('Phase_',num2str(center_freq));
    phaseandfreq = strcat(phaseandfreq, '_');
    mat_name = strcat(phaseandfreq, mat_name);
    mat_name = strcat(mat_name,'.mat');
    disp(['Saving file:' mat_name ' in folder:' mat_pathstr]);
    save(mat_name);
    %save(mat_name, 'corrmat', 'corr_mat_name');  
    disp('DONE!')
     
    figure
    disp([ 'Drawing color matrix for phase sync....' ])
    %subplot(221)
    %t=({['ISPC, patient' num2str(x)];['text_B=' num2str(y)];['text_C=' num2str(y)]});
    clims = [0 1];
    imagesc(phase_synchronization_matrix,clims);
    colorbar;
    t=({['Phase-based synchronicity, Patient:' patient ' ,Frequency=' num2str(center_freq)];['Condition:' condition]});
    title(t);
    xlabel('Channel') % x-axis label
    ylabel('Channel') % y-axis label});
    title(t);
    xlabel({[ ' Mean polar abs=' num2str(polar_abs) ', Mean polar angle =' num2str(polar_angle)]}) % x-axis label
    ylabel('Channel') % y-axis label
    disp({[ ' Mean polar abs=' num2str(polar_abs) ', Mean polar angle =' num2str(polar_angle)]})
    %draw_chart_polar (patient, condition)
end
end

