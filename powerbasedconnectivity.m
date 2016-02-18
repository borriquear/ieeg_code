%% This function does powerbased connectivity computation -ch27-

function [] = powerbasedconnectivity()
%% powerbasedconnectivity  Calculate r_pearson_corrmat and r_spearman_corrmat, savein mat files if all_pair= 1 AND chantouse = EEG.nbchan
% INPUT: centefreq
% shows the distribution of powerdata per channel, both raw and logartihm
% scale to make it look more normal. We do the SK test of normality
% can be used for single channels or for all channels, iof single channel
% draw disribution, f all channels only if chanel power time series has
% normal distribution
% 1. Load data
%    Call to initialize_EEG_variables() and initialize_wavelet(EEG_study)
% 2. Calculate if the channels have a normal distribution, if so, display the power time series and the log power
    % We check to display the chart only for those channels that have no normal
    % distribution, or when we are interested in only one channel, initchan = chantouse
% 3. Calculate the correlation matrix between pairs or only one pair(all_pair)
    % all_pair= 1, default,  to calculate corretaltion for all pairs, 
    % r_pearson_corrmat and r_spearman_corrmat
% 4. Call to timefrequencypower(...) to calculate Time frequency power, for 2 channels conditional to another
% 
%% 1. Load epoch and Initialize data
%load sampleEEGdata
%disp('Closing open windows....')
%close all
disp('Loading data....')
[myfullname, EEG_study, channel_labels] = initialize_EEG_variables();
EEG = EEG_study;
[half_of_wavelet_size, n_wavelet, n_data, n_convolution, wavelet_cycles, time] = initialize_wavelet(EEG_study);
% draw_all_single_correlations = 0 do not draw all correlation chart,
% justmatrix, set to 0 when calculating the entire matrix
draw_all_single_correlations = 0;
times2plot =  dsearchn(EEG.times',[ EEG.times(1) EEG.times(end)]');
trial2plot = EEG.trials;
[pathstr, subjandcond, matext] = fileparts(myfullname);
%To plot only one channel: initchan = chantouse
%To plot from channel 12 to 23, initchan = 12, chantouse = 23;
initchan =  2;
%chantouse = 11;
chantouse = EEG.nbchan;
%chantouse = 3;
disp('Data Loaded OK!');
center_freq_bands = [2.5 5 10 20 40] %delta, theta, alpha,beta, gamma
%center_freq_bands = [2.5]
disp(['Calculating Power-based connectibity for ' num2str(size(center_freq_bands,2)) ' frequency bands '] )
%% 2..4
for freq_index=1:size(center_freq_bands,2)
    centerfreq = center_freq_bands(freq_index); 
    %% 2. Calculate if the channels have a normal distribution, if so, display the power time series and the log power
    % Show the chart only for one channel
    disp('Calculating if there is some channel with power time series Nornally distributed...')
    %disp(['Calculating for channel:', num2str(initchan) ' to ', num2str(chantouse)])
    %number of channels with power values normally distributed
    counternormalc = 0;
    
    for chani=initchan:chantouse
        %setup wavelet convolution and outputs
        sensor2use = channel_labels(chani);      
        % FFT of data (note: this doesn't change on frequency iteration)
        fft_data = fft(reshape(EEG.data(strcmpi(sensor2use,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);
        
        % create wavelet and run convolution
        fft_wavelet            = fft(exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq))^2)),n_convolution);
        convolution_result_fft = ifft(fft_wavelet.*fft_data,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq));
        convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
        convolution_result_fft = abs(reshape(convolution_result_fft,EEG.pnts,EEG.trials)).^2;
        
        % trim edges so the distribution is not driven by edge artifact outliers
        % (note: here we just use visual inspection to remove edges)
        convolution_result_fft = convolution_result_fft(100:end-100,:);
        % test for normal distribution, if you have the stats toolbox
        if exist('kstest','file')
            [h,p1] = kstest(convolution_result_fft(:));
            [h,p2] = kstest(log10(convolution_result_fft(:)));
            [h,p3] = kstest(randn(numel(convolution_result_fft),1));
            disp([ 'Channel:' sensor2use 'KS test for norm power: '        num2str(p1) ' (>.05 means normal distribution) ' ])
            disp([ 'Channel:' sensor2use 'KS test for log norm power: ' num2str(p2) ' (>.05 means normal distribution) ' ])
            %disp([ 'KS test for normality of random data: '  num2str(p3) ' (>.05 means normal distribution) ' ])
        end
        % display the chart only for those channels that have no normal
        % distribution, or when we are interested in only one channel, initchan =
        % chantouse
        if (p1 + p2 > 0.05) || (initchan == chantouse)
            % plot distirbution of power data
            %counter of channel with normal distribtion
            counternormalc = counternormalc  + 1;
            figure;
            %subplot(chantouse-1,2,(chani-1)*2-1)
            subplot(1,2,1)
            hist(convolution_result_fft(:),500)
            % str=sprintf('Distribution of power values = %s ', num2str(sensor2use));
            title(['Distribution of power for channel: ' sensor2use ])
            %title(str)
            %title('Distribution of power values' sensor2use )
            axis square
            %set(gca,'xlim',[EEG.times(1) EEG.times(end)],'ylim',[min(real_data(:)) max(real_data(:))])
            xlabel(['Power (\muV)) channel:  ' sensor2use ])
            ylabel('Count')
            %subplot(chantouse-1,2,(chani-1)*2)*
            subplot(1,2,2)
            hist(log10(convolution_result_fft(:)),500)
            title(['Distribution of log_1_0power for channel:', sensor2use])
            xlabel(['log10 Power for channel: ' sensor2use ])
            ylabel('Count')
            axis square
        end
    end
    
    %% 3. Calculate the correlation matrix between pairs or only one pair
    % all_pair = 1 corr matrix  for all combinations between initchan and chantouse, all_pair = 0 only between two channels
    all_pair= 1;
    if all_pair == 0
        % Calculate the correlation of only one pair
        sensor1 = 'LHD2';
        sensor2 = 'LHD1';
        disp(['Estimation of correlation coefficient of pair ' sensor1 '--' sensor2])
        %call compute_convolution
        [convolution_result_fft1, convolution_result_fft2,r_pearson, r_spearman,] = compute_convolution_onepair(sensor1, sensor2, EEG, half_of_wavelet_size, n_wavelet, n_data, n_convolution, wavelet_cycles, time, centerfreq, trial2plot,times2plot, draw_all_single_correlations);
        disp(['Estimation of correlation coefficient r_pearson=' num2str(r_pearson) ' , r_spearman=' num2str(r_spearman) ])
    else
        %initialize matrix of correlation coefficients
        %number of rows (=cols) for which is actually calculated the correlation the correlation matrix
        nbofrows = chantouse - initchan + 1;
        r_pearson_corrmat = zeros(EEG.nbchan-1,EEG.nbchan-1);
        r_spearman_corrmat = zeros(EEG.nbchan-1,EEG.nbchan-1);
        disp(['Estimation of correlation coefficient of ALL pair between' channel_labels(initchan) '--' channel_labels(chantouse)])
        for chani=initchan:chantouse
            for chanj=chani +1:chantouse
                sensor1 = channel_labels(chani);
                sensor2 = channel_labels(chanj);
                [convolution_result_fft1, convolution_result_fft2,r_pearson, r_spearman] = compute_convolution_onepair(sensor1, sensor2, EEG, half_of_wavelet_size, n_wavelet, n_data, n_convolution, wavelet_cycles, time, centerfreq, trial2plot,times2plot,draw_all_single_correlations);
                disp(['Estimation of correlation coefficient of' sensor1 '--' sensor2 ' r_pearson=' num2str(r_pearson) ' r_spearman=' num2str(r_spearman) ' DONE!'])
                r_pearson_corrmat(chani-1,chanj-1) = r_pearson;
                r_pearson_corrmat(chanj-1,chani-1) = r_pearson_corrmat(chani-1,chanj-1);
                r_spearman_corrmat(chani-1,chanj-1) = r_spearman;
                r_spearman_corrmat(chanj-1,chani-1) = r_spearman_corrmat(chani-1,chanj-1);
            end
        end
        %calculate correlation matrix
        %save r_spearman_corrmat as mat file and draw colour matrix
        typeofcorr = 'Pearson';
        save_and_display_corr_mat(r_pearson_corrmat, initchan, nbofrows, channel_labels,centerfreq, wavelet_cycles, subjandcond, typeofcorr,myfullname);
        typeofcorr = 'Spearman';
        save_and_display_corr_mat(r_spearman_corrmat, initchan, nbofrows, channel_labels,centerfreq, wavelet_cycles, subjandcond, typeofcorr,myfullname);
    end
    disp(['DONE for frequency =' num2str(centerfreq)])
end
%% 4. Time frequency power, for 2 channels and 2 channels conditional to another
disp('Calling to time-freq power for two channels....')
%res_func = timefrequencypower(subjandcond, EEG, centerfreq, time, half_of_wavelet_size, n_wavelet, n_data, n_convolution, wavelet_cycles);

end

function [convolution_result_fft1, convolution_result_fft2, r_pearson, r_spearman] = compute_convolution_onepair(sensor1, sensor2, EEG, half_of_wavelet_size, n_wavelet, n_data, n_convolution, wavelet_cycles, time,centerfreq, trial2plot,times2plot, draw_all_single_correlations)
%% Figure 27.4
fft_data1 = fft(reshape(EEG.data(strcmpi(sensor1,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);
fft_data2 = fft(reshape(EEG.data(strcmpi(sensor2,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);
% create wavelet and run convolution
fft_wavelet             = fft(exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq))^2)),n_convolution);
convolution_result_fft  = ifft(fft_wavelet.*fft_data1,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq));
convolution_result_fft  = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
convolution_result_fft1 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);

fft_wavelet             = fft(exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq))^2)),n_convolution);
convolution_result_fft  = ifft(fft_wavelet.*fft_data2,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq));
convolution_result_fft  = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
convolution_result_fft2 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);

% keep only requested time regions
convolution_result_fft1 = convolution_result_fft1(times2plot(1):times2plot(2),:);
convolution_result_fft2 = convolution_result_fft2(times2plot(1):times2plot(2),:);
r_pearson = corr(abs(convolution_result_fft1(:,trial2plot)).^2,abs(convolution_result_fft2(:,trial2plot)).^2,'type','p');
r_spearman = corr(abs(convolution_result_fft1(:,trial2plot)).^2,abs(convolution_result_fft2(:,trial2plot)).^2,'type','s');

%call to the function to draw the charts for the correlation between
%sensor1 and sensor2
drawcorrelationchart_onepair(r_pearson, r_spearman, EEG,convolution_result_fft1, convolution_result_fft2, sensor1, sensor2, times2plot, trial2plot, centerfreq, draw_all_single_correlations);
end

function [] = drawcorrelationchart_onepair(r_pearson, r_spearman,EEG,convolution_result_fft1, convolution_result_fft2, sensor1, sensor2, times2plot, trial2plot, centerfreq, draw_all_single_correlations)
%% Draw correlation between the given pair of channels
% if draw_all_single_correlations = 1 it actually  draws the chart
% otherwise it does not
%global draw_all_single_correlations
% Draw all the correlation charts for each combination of pair of channels,
% not a good idea, NxN charts
if draw_all_single_correlations == 1
    figure
    subplot(211)
    plot(EEG.times(times2plot(1):times2plot(2)),abs(convolution_result_fft1(:,trial2plot)).^2)
    hold on
    plot(EEG.times(times2plot(1):times2plot(2)),abs(convolution_result_fft2(:,trial2plot)).^2,'r')
    xlabel('Time (ms)')
    set(gca,'xlim',EEG.times(times2plot))
    %legend({sensor1;sensor2})
    
    subplot(223)
    plot(abs(convolution_result_fft1(:,trial2plot)).^2,abs(convolution_result_fft2(:,trial2plot)).^2,'.')
    title('Person: Power relationship')
    xlabel([ sensor1 ' ' num2str(centerfreq) 'Hz power' ])
    ylabel([ sensor2 'Hz power' ])
    %r_pearson=corr(abs(convolution_result_fft1(:,trial2plot)).^2,abs(convolution_result_fft2(:,trial2plot)).^2,'type','p');
    legend([ 'Pearson R = ' num2str(r_pearson) ' at Freq=' num2str(centerfreq) ]);
    
    subplot(224)
    plot(tiedrank(abs(convolution_result_fft1(:,trial2plot)).^2),tiedrank(abs(convolution_result_fft2(:,trial2plot)).^2),'.')
    title('Spearman or Rank-power relationship')
    xlabel([ sensor1 ' ' num2str(centerfreq) 'Hz rank-power' ])
    ylabel([ sensor2  'Hz rank-power' ])
    %r_spearman=corr(abs(convolution_result_fft1(:,trial2plot)).^2,abs(convolution_result_fft2(:,trial2plot)).^2,'type','s');
    legend([ 'Spearman Rho = ' num2str(r_spearman)  ' at Freq=' num2str(centerfreq)]);
    set(gca,'ylim',get(gca,'xlim'))
end
end

function [ret_val] =  save_and_display_corr_mat(corrmat, initialrow, nbofrows, channel_labels, centerfreq, wavelet_cycles, subjandcond, typeofcorr, myfullname)
%% Saves as a mat file and draws correlation matrix with contour using imagesc
% from the entire corrmat it cuts a section given from initial to chantouse
% and draws that sub matrix

ret_val = 0;
sqcorrmatsize = size(corrmat);
totalnbrows = sqcorrmatsize(1);
if issym(corrmat) == 0
    msgerror = 'ERROR The correlation matrix is not symmetric';
    error(msgerror);
    ret_val = -1;
else
    %Save object  corrmat as mat file
    [mat_pathstr,mat_name,mat_ext] = fileparts(myfullname);
    %typeofcorr =strcat(typeofcorr, '_');
    ftypeofcorr =strcat(typeofcorr, '_Fq_');
    %freq = sprintf('%d_', centerfreq);
    ftypeofcorr =strcat(ftypeofcorr,num2str(centerfreq));
    ftypeofcorr = strcat(ftypeofcorr,'_');
    corr_mat_name = strcat(ftypeofcorr,mat_name);
    %create directory if does not exist 
    cd(mat_pathstr);
    if ~exist(mat_name, 'dir')
        mkdir(mat_name);
    end
    cd(mat_name);
    mat_pathstr = strcat(mat_pathstr, '\');
    target_dir_mat = strcat(mat_pathstr, mat_name);
    corr_mat_name = strcat(corr_mat_name, '.mat');
    corr_mat_name = fullfile(target_dir_mat,corr_mat_name);
    disp(['Saving correlation matrix in ...' corr_mat_name]);
    save(corr_mat_name, 'corrmat', 'corr_mat_name');  
    disp('DONE!')
    disp('Displaying the corr matrix...')
    figure;
    %Build corrmat for the channels that is being actually calculated
    %initialize matrix nbofrows x nbofrows to zeros
    corrmat_reduced = zeros(nbofrows, nbofrows);
    finalrow = nbofrows + initialrow -1;
    auxi = 1;
    auxj = 1;
    strofchannels = {};
    for i=initialrow:finalrow
        strofchannels = [strofchannels, channel_labels(i)];
    end
    for i=initialrow-1:finalrow-1
        for j =initialrow-1:finalrow-1
            corrmat_reduced(auxi,auxj) = corrmat(i,j);
            corrmat_reduced(auxj,auxi) = corrmat_reduced(auxi,auxj); 
            auxj = auxj +1;   
        end
        auxj = 1;
        auxi = auxi + 1;
    end
    %strofchannels = {'sam'; 'alan'; 'ellie'};

    imagesc(corrmat_reduced);
    ax = gca;
    ax.XTickLabelRotation=45;
    set(gca, 'XTickLabel',strofchannels, 'YTickLabel',strofchannels,'XTick',1:numel(strofchannels),'YTick',1:numel(strofchannels))   
    
    %typeofcorr = 'Spearman';
    msgtitle = sprintf('Correlation Coefficient between pairs: %s ', typeofcorr);
    title(msgtitle);
    texstring = texlabel(subjandcond,'literal');
    %texstring = strsplit(texstring, '\_');
    %xlabel(['Cond: ', texstring(4), ' ,Patient=', texstring(5), 'Date session=' texstring(6),texstring(7), ' ,at Frequency=' num2str(centerfreq), ' ,Wavelet #cycles=' num2str(wavelet_cycles) ])
    xlabel([' ', texstring, ' ,at Frequency=' num2str(centerfreq), ' ,Wavelet #cycles=' num2str(wavelet_cycles) ])
    colormap('jet');
    colorbar;
    
end
end

%% Figure 27.7
%load sampleEEGdata 
function [res_func] = timefrequencypower(subjandcond, EEG, centerfreq, time, half_of_wavelet_size, n_wavelet, n_data, n_convolution, wavelet_cycles)
%% time frequency power calculation between two channels and a cotrol channel
disp('Calculating  time frequency map of correlation for two channels, accross frequencies')

seed_chan     = 'LPT1';
target_chan   = 'LPT5';
control_chan  = 'RPT1';
disp(sprintf('Seed chan = %s', seed_chan))
disp(sprintf('Target chan = %s', target_chan))
disp(sprintf('Control chan = %s', control_chan))
clim = [0 1];

% wavelet parameters
min_freq = 2;
max_freq = 40;
num_frex = 15;
disp(sprintf('Range of frequencies: min-max-total = %d , %d , %d', min_freq, max_freq, num_frex))
res_func = 1;

% downsampled times
%times2save = -200:50:800;
% times2save = EEG.times; % uncomment this line for figure 27.8
times2save = 1000:1:EEG.times(end);
%times2save = EEG.times; 

% other wavelet parameters
frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
%time = -1:1/EEG.srate:1;
%half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters
%n_wavelet     = length(time);
%n_data        = EEG.pnts*EEG.trials;
%n_convolution = n_wavelet+n_data-1;
%wavelet_cycles= 4.5;

times2saveidx = dsearchn(EEG.times',times2save');

% FFT of data (note: this doesn't change on frequency iteration)
fft_data_seed = fft(reshape(EEG.data(strcmpi(seed_chan,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);
fft_data_trgt = fft(reshape(EEG.data(strcmpi(target_chan,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);
fft_data_ctrl = fft(reshape(EEG.data(strcmpi(control_chan,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);


% initialize output time-frequency data
tf_corrdata = zeros(length(frequencies),length(times2save),2);
disp('Calculating fft and correlation coefficient accross all frequencies')
for fi=1:length(frequencies)
    
    % create wavelet and get its FFT
    fft_wavelet = fft(exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi),n_convolution);
    
    
    % convolution for all three sites (save only power)
    convolution_result_fft = ifft(fft_wavelet.*fft_data_seed,n_convolution) * sqrt(wavelet_cycles /(2*pi*frequencies(fi)));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    conv_result_seed = abs(reshape(convolution_result_fft,EEG.pnts,EEG.trials)).^2;
    
    convolution_result_fft = ifft(fft_wavelet.*fft_data_trgt,n_convolution) * sqrt(wavelet_cycles /(2*pi*frequencies(fi)));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    conv_result_trgt = abs(reshape(convolution_result_fft,EEG.pnts,EEG.trials)).^2;
    
    convolution_result_fft = ifft(fft_wavelet.*fft_data_ctrl,n_convolution) * sqrt(wavelet_cycles /(2*pi*frequencies(fi)));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    conv_result_ctrl = abs(reshape(convolution_result_fft,EEG.pnts,EEG.trials)).^2;
    
    % downsample and rank transform all data
    conv_result_seed = tiedrank(conv_result_seed(times2saveidx,:)')';
    conv_result_trgt = tiedrank(conv_result_trgt(times2saveidx,:)')';
    conv_result_ctrl = tiedrank(conv_result_ctrl(times2saveidx,:)')';
    
    for ti=1:length(times2save)
        denomtrials = EEG.trials*(EEG.trials^2-1);
        if denomtrials == 0
            denomtrials = 1;
        end    
        % compute bivariate correlations
        r_st = 1-6*sum((conv_result_seed(ti,:)-conv_result_trgt(ti,:)).^2)/denomtrials;
        r_sc = 1-6*sum((conv_result_seed(ti,:)-conv_result_ctrl(ti,:)).^2)/denomtrials;
        r_tc = 1-6*sum((conv_result_ctrl(ti,:)-conv_result_trgt(ti,:)).^2)/denomtrials;
        
        % bivariate correlation for comparison
        tf_corrdata(fi,ti,1) = r_st;

        % compute partial correlation and store in results matrix
        tf_corrdata(fi,ti,2) = (r_st-r_sc*r_tc) / ( sqrt(1-r_sc^2)*sqrt(1-r_tc^2) );
    end
    disp(sprintf('...calculated corr. map for frequency =%s', num2str(frequencies(fi))))
end
disp('Displaying time-frequency map of correlation between 2 channels')
figure
for i=1:2
    subplot(1,2,i)
    contourf(times2save,frequencies,squeeze(tf_corrdata(:,:,i)),40,'linecolor','none')
    set(gca,'clim',clim,'xlim',[times2save(1) times2save(end)],'yscale','log','ytick',logspace(log10(frequencies(1)),log10(frequencies(end)),6),'yticklabel',round(logspace(log10(frequencies(1)),log10(frequencies(end)),6)*10)/10)
    axis square
    if i==1
        title([ 'Spearman Correlation ' seed_chan ',' target_chan ])
    else
        title([ 'Spearman Partial correlation ' seed_chan ',' target_chan ' | ' control_chan ])
    end
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    
    %texstring = strsplit(texstring, '\_');
    %xlabel(['Cond: ', texstring(4), ' ,Patient=', texstring(5), 'Date session=' texstring(6),texstring(7), ' ,at Frequency=' num2str(centerfreq), ' ,Wavelet #cycles=' num2str(wavelet_cycles) ])
end
texstring = texlabel(subjandcond,'literal');
legend(['', texstring, ' ,Wavelet #cycles=' num2str(wavelet_cycles) ])
res_func = 1;
end
