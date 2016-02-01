%% This function does powerbased connectivity computation -ch27-

function [] = powerbasedconnectivity()
%% Figure 27.2
% shows the distribution of powerdata per channel, both raw and logartihm
% scale to make it look more normal. We do the SK test of normality
% can be used for single channels or for all channels, iof single channel
% draw disribution, f all channels only if chanel power time series has
% normal distribution
%load sampleEEGdata
%% Load epcoh and Initialize data
%
disp('Closing open windows....')
close all
disp('DONE! Closing open windows....')
disp('Loading data....')
[myfullname, EEG_study, channel_labels, time] = initialize_EEG_variables();
[half_of_wavelet_size, n_wavelet, n_data n_convolution, wavelet_cycles] = initialize_wavelet(time, EEG_study);
EEG = EEG_study;
draw_all_single_correlations = 0;
times2plot =  dsearchn(EEG.times',[ EEG.times(1) EEG.times(end)]');
trial2plot = EEG.trials;
centerfreq = 10;
[pathstr,subjandcond,matext] = fileparts(myfullname);
% Which channels we want to calculate the correlation coefficient
%To plot only one channel: initchan = chantouse
%To plot from channel 12 to 23, initchan = 12, chantouse = 23;
initchan =  2;
chantouse = 11;
chantouse = EEG.nbchan;
chantouse = 3;
disp('Data Loaded OK!');

%% Calculate if the channels have a normal distribution, of so, display the power time series and the log
% Show the chart only for one channel 
disp('Calculating if there is some channel with power time series Nornally distributed...')
disp(['Calculating for channel:', num2str(initchan) ' to ', num2str(chantouse)])
%number of channels with power values normally distributed
counternormalc = 0;
for chani=initchan:chantouse 
    %setup wavelet convolution and outputs
    sensor2use = channel_labels(chani);
    time = -1:1/EEG.srate:1;
    half_of_wavelet_size = (length(time)-1)/2;
    % FFT parameters
    n_wavelet     = length(time);
    n_data        = EEG.pnts*EEG.trials;
    n_convolution = n_wavelet+n_data-1;
    wavelet_cycles= 4.5;
    
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
%% Calculate the correlation matrix between pairs or only one pair
% all_pair = 1 corr matrix  for all combinations between initchan and chantouse, all_pair = 0 only between two channels
all_pair= 1;
if all_pair == 0 
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
    save_and_display_corr_mat(r_spearman_corrmat, initchan, nbofrows, channel_labels,centerfreq, wavelet_cycles, subjandcond, typeofcorr);
    typeofcorr = 'Spearman';
    save_and_display_corr_mat(r_spearman_corrmat, initchan, nbofrows, channel_labels,centerfreq, wavelet_cycles, subjandcond, typeofcorr);
    
end
%% Figure 27.5
return
% Compute how many time points are in one cycle, and limit xcov to this lag
nlags = round(EEG.srate/centerfreq);
%% 

% note that data are first tiedrank'ed, which results in a Spearman rho
% instead of a Pearson r. 
[corrvals,corrlags] = xcov(tiedrank(abs(convolution_result_fft1(:,trial2plot)).^2),tiedrank(abs(convolution_result_fft2(:,trial2plot)).^2),nlags,'coeff');

% convert correlation lags from indices to time in ms
corrlags = corrlags * 1000/EEG.srate; 

figure
plot(corrlags,corrvals,'-o','markerface','w')
hold on
plot([0 0],get(gca,'ylim'))
xlabel([ sensor1 ' leads --- Time lag in ms --- ' sensor2 ' leads' ])
ylabel('Correlation coefficient')

%% Figure 27.6
% this is over trials, willnot work.

load sampleEEGdata 
sensor1 = 'poz';
sensor2 = 'fz';
%sensor1 = 'RPT1';
%sensor2 = 'RPT2';

timewin1 = [ -300 -100 ]; % in ms relative to stim onset
timewin2 = [  200  400 ];

centerfreq1 =  6; % in Hz
centerfreq2 =  6;


% convert time from ms to index
timeidx1 = zeros(size(timewin1));
timeidx2 = zeros(size(timewin2));
for i=1:2
    [junk,timeidx1(i)] = min(abs(EEG.times-timewin1(i)));
    [junk,timeidx2(i)] = min(abs(EEG.times-timewin2(i)));
end

% setup wavelet convolution and outputs
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
wavelet_cycles= 4.5;

% FFT of data (note: this doesn't change on frequency iteration)
fft_data1 = fft(reshape(EEG.data(strcmpi(sensor1,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);
fft_data2 = fft(reshape(EEG.data(strcmpi(sensor2,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);

% initialize output time-frequency data
corrdata = zeros(EEG.trials,2);

% create wavelet and run convolution
fft_wavelet            = fft(exp(2*1i*pi*centerfreq1.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq1))^2)),n_convolution);
convolution_result_fft = ifft(fft_wavelet.*fft_data1,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq1));
convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
analyticsignal1        = abs(convolution_result_fft).^2;

fft_wavelet            = fft(exp(2*1i*pi*centerfreq2.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq2))^2)),n_convolution);
convolution_result_fft = ifft(fft_wavelet.*fft_data2,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq2));
convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
analyticsignal2        = abs(convolution_result_fft).^2;


% Panel A: correlation in a specified window
tfwindowdata1 = mean(analyticsignal1(timeidx1(1):timeidx1(2),:),1);
tfwindowdata2 = mean(analyticsignal2(timeidx2(1):timeidx2(2),:),1);
figure
subplot(121)
plot(tfwindowdata1,tfwindowdata2,'.')
axis square
title([ 'TF window correlation, r_p=' num2str(corr(tfwindowdata1',tfwindowdata2','type','p')) ])
xlabel([ sensor1 ': ' num2str(timewin1(1)) '-' num2str(timewin1(2)) '; ' num2str(centerfreq1) ' Hz' ])
ylabel([ sensor2 ': ' num2str(timewin2(1)) '-' num2str(timewin2(2)) '; ' num2str(centerfreq2) ' Hz' ])


% also plot rank-transformed data
subplot(122)
plot(tiedrank(tfwindowdata1),tiedrank(tfwindowdata2),'.')
axis square
xlabel([ sensor1 ': ' num2str(timewin1(1)) '-' num2str(timewin1(2)) '; ' num2str(centerfreq1) ' Hz' ])
ylabel([ sensor2 ': ' num2str(timewin2(1)) '-' num2str(timewin2(2)) '; ' num2str(centerfreq2) ' Hz' ])
title([ 'TF window correlation, r_p=' num2str(corr(tfwindowdata1',tfwindowdata2','type','s')) ])


% panel B: correlation over time
corr_ts = zeros(size(EEG.times));
for ti=1:EEG.pnts
    corr_ts(ti) = corr(analyticsignal1(ti,:)',analyticsignal2(ti,:)','type','s');
end

figure
plot(EEG.times,corr_ts)
set(gca,'xlim',[-200 1200])
xlabel('Time (ms)'), ylabel('Spearman''s rho')



% Panel C: exploratory time-frequency power correlations
times2save = -200:25:1200;
frex = logspace(log10(2),log10(40),20);


times2save_idx = zeros(size(times2save));
for i=1:length(times2save)
    [junk,times2save_idx(i)] = min(abs(EEG.times-times2save(i)));
end

% rank-transforming the data can happen outside the frequency loop
seeddata_rank = tiedrank(tfwindowdata2);

% initialize output correlation matrix
expl_corrs = zeros(length(frex),length(times2save));

for fi=1:length(frex)
    
    % get power (via wavelet convolution) from signal1
    fft_wavelet            = fft(exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frex(fi)))^2)),n_convolution);
    convolution_result_fft = ifft(fft_wavelet.*fft_data1,n_convolution) * sqrt(wavelet_cycles /(2*pi*frex(fi)));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
    analyticsignal1        = abs(convolution_result_fft).^2;
    
    for ti=1:length(times2save)
        expl_corrs(fi,ti) = 1-6*sum((seeddata_rank-tiedrank(analyticsignal1(times2save_idx(ti),:))).^2)/(EEG.trials*(EEG.trials^2-1));
    end
end

figure
contourf(times2save,frex,expl_corrs,40,'linecolor','none')
set(gca,'clim',[-.4 .4],'yscale','log','ytick',round(logspace(log10(frex(1)),log10(frex(end)),8)))
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'Correlation over trials from seed ' sensor2 ', ' num2str(centerfreq2) ' Hz and ' num2str(timewin2(1)) '-' num2str(timewin2(2)) ' ms' ])
colorbar

%% Figure 27.7
%load sampleEEGdata 
[myfullname, EEG_study, channel_labels, initchan, chantouse, time] = initialize_EEG_variables();
[half_of_wavelet_size, n_wavelet, n_data n_convolution, wavelet_cycles] = initialize_wavelet(time, EEG_study);
EEG = EEG_study;
seed_chan     = 'fz';
target_chan   = 'f6';
control_chan  = 'f1';
seed_chan     = 'LPT1';
target_chan   = 'LPT5';
control_chan  = 'RPT1';

clim = [0 1];

% wavelet parameters
min_freq = 2;
max_freq = 40;
num_frex = 15;

% downsampled times
times2save = -200:50:800;
% times2save = EEG.times; % uncomment this line for figure 27.8
times2save = 0:50:10000;
times2save = EEG.times; 

% other wavelet parameters
frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
wavelet_cycles= 4.5;

times2saveidx = dsearchn(EEG.times',times2save');

% FFT of data (note: this doesn't change on frequency iteration)
fft_data_seed = fft(reshape(EEG.data(strcmpi(seed_chan,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);
fft_data_trgt = fft(reshape(EEG.data(strcmpi(target_chan,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);
fft_data_ctrl = fft(reshape(EEG.data(strcmpi(control_chan,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);


% initialize output time-frequency data
tf_corrdata = zeros(length(frequencies),length(times2save),2);

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
end

% plot
figure
for i=1:2
    subplot(1,2,i)
    contourf(times2save,frequencies,squeeze(tf_corrdata(:,:,i)),40,'linecolor','none')
    set(gca,'clim',clim,'xlim',[-200 800],'yscale','log','ytick',logspace(log10(frequencies(1)),log10(frequencies(end)),6),'yticklabel',round(logspace(log10(frequencies(1)),log10(frequencies(end)),6)*10)/10)
    axis square
    if i==1
        title([ 'Correlation between ' seed_chan ' and ' target_chan ])
    else
        title([ 'Partial correlation between ' seed_chan ' and ' target_chan ])
    end
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
end

%% Figure 27.8

% Re-run the code for the previous figure but comment out the 
% following line towards the top:
% times2save = EEG.times; % uncomment this line for figure 27.8
% Then run this section of code.


ds_timesidx = dsearchn(EEG.times',(-200:50:800)'); % ds = downsampled
[~,lofreq]  = min(abs(frequencies-4.7));
[~,hifreq]  = min(abs(frequencies-32));

figure

subplot(221)
contourf(times2save,frequencies,squeeze(tf_corrdata(:,:,2)),40,'linecolor','none')
hold on
plot(get(gca,'xlim'),frequencies([lofreq lofreq]),'k--')
plot(get(gca,'xlim'),frequencies([hifreq hifreq]),'k--')
set(gca,'clim',clim,'xlim',[-200 800],'yscale','log','ytick',logspace(log10(frequencies(1)),log10(frequencies(end)),6),'yticklabel',round(logspace(log10(frequencies(1)),log10(frequencies(end)),6)*10)/10)
title('Original (256 Hz)')


subplot(222)
contourf(times2save(ds_timesidx),frequencies,squeeze(tf_corrdata(:,ds_timesidx,2)),40,'linecolor','none')
hold on
plot(get(gca,'xlim'),frequencies([lofreq lofreq]),'k--')
plot(get(gca,'xlim'),frequencies([hifreq hifreq]),'k--')
set(gca,'clim',clim,'xlim',[-200 800],'yscale','log','ytick',logspace(log10(frequencies(1)),log10(frequencies(end)),6),'yticklabel',round(logspace(log10(frequencies(1)),log10(frequencies(end)),6)*10)/10)
title('Down-sampled (20 Hz)')


subplot(223)
plot(EEG.times,squeeze(tf_corrdata(lofreq,:,2)))
hold on
plot(EEG.times(ds_timesidx),squeeze(tf_corrdata(lofreq,ds_timesidx,2)),'ro-','markerface','w')
title('Effect of downsampling on low-frequency activity')
set(gca,'xlim',[-200 800],'ylim',[.25 .65])


subplot(224)
plot(EEG.times,squeeze(tf_corrdata(hifreq,:,2)))
hold on
plot(EEG.times(ds_timesidx),squeeze(tf_corrdata(hifreq,ds_timesidx,2)),'ro-','markerface','w')
title('Effect of downsampling on high-frequency activity')
set(gca,'xlim',[-200 800],'ylim',[-.1 .6])
legend({'Original (256 Hz)';'Down-sampled (20 Hz)'})

%% Figure 27.9

% note: this cell takes a while to run, particularly on slow computers!

n = 1000;
ncorrs = 100000;

t=[0 0 0];
for i=1:ncorrs
    
    % create random variables
    a = rand(2,n);
    
    tic
    % Matlab corr function
    c2 = corr(a','type','s');
    t(1) = t(1) + toc;
    
    
    tic
    % self-written Spearman correlation (must first rank-transform)
    a1 = tiedrank(a')'; % tiedrank accepts matrix input, but make sure the matrix is in the correct orientation!!
    c1 = 1-6*sum((a1(1,:)-a1(2,:)).^2)/(n*(n^2-1));
    t(2) = t(2) + toc;
    
    
    tic
    % ordinary least squares
    % Note: Uncommenting the following line will normalize the data to give 
    % you a correlation coefficient. If you don't need the correlation coefficient 
    % (and instead can use unstandardized regression coefficients), leave this
    % line commented out for a ten-fold increase in speed.
    %a = bsxfun(@rdivide,bsxfun(@minus,a,mean(a,2)),std(a,[],2));
    c3 = (a(1,:)*a(1,:)')\a(1,:)*a(2,:)';
    t(3) = t(3) + toc;
end

figure
bar(t)
set(gca,'xticklabel',{'corr function';'manual';'ols'},'xlim',[.5 3.5])
ylabel([ 'Time for ' num2str(ncorrs) ' iterations (s)' ])

%% end.
%%% Not Useful charts  %%%%%
%% 27.0
%% initial example to calculate correlation between 2 random variables using different functions
%NOP
% corr1 is Matlab corrcoef, coor2 calculating the formula ourselves
% covariance scaled by vaiance, corr3 with matrixes
a = randn(1,100);
b = randn(1,100);

corr1 = corrcoef(a,b);

a1=a-mean(a);
b1=b-mean(b);
corr2 = (a1*b1')/sqrt( (a1*a1')*(b1*b1') );

c = [a1; b1];
covmat = c*c';

% notice the following:
covmat(1,1) == a1*a1'
covmat(2,2) == b1*b1'
covmat(2,1) == a1*b1' 
% actually, some of these might not be exactly equal due to very small computer rounding errors. 
% try this instead:
(covmat(2,1)-a1*b1') < 0.0000000000001

corr3 = covmat(1,2)/sqrt(covmat(1)*covmat(end));

fprintf([ '\nMatlab corrcoef function: ' num2str(corr1(1,2)) '\n' ...
            'covariance scaled by variances: ' num2str(corr2) '\n' ...
            'covariance computed as matrix: ' num2str(corr3) '\n\n' ])

%% Figure 27.1
%NOP
% to show the difference between pearson and spearman
% when data are normally distribted, P and S are almost identical, but if
% there are outliers or the distribution is nonnormal then P inflates or
% deflates the correlation
anscombe = [
  % series 1    series 2    series 3     series 4
    10  8.04    10  9.14    10  7.46     8  6.58;
     8  6.95     8  8.14     8  6.77     8  5.76;
    13  7.58    13  8.76    13 12.74     8  7.71;
     9  8.81     9  8.77     9  7.11     8  8.84;
    11  8.33    11  9.26    11  7.81     8  8.47;
    14  9.96    14  8.10    14  8.84     8  7.04;
     6  7.24     6  6.13     6  6.08     8  5.25;
     4  4.26     4  3.10     4  5.39     8  5.56;
    12 10.84    12  9.13    12  8.15     8  7.91;
     7  4.82     7  7.26     7  6.42     8  6.89;
     5  5.68     5  4.74     5  5.73    19 12.50;
    ];


% plot and compute correlations
figure
for i=1:4
    subplot(2,2,i)
    %accurate, inflated, deflated , inflated
    plot(anscombe(:,(i-1)*2+1),anscombe(:,(i-1)*2+2),'.')
    lsline
    corr_p = corr(anscombe(:,(i-1)*2+1),anscombe(:,(i-1)*2+2),'type','p');
    corr_s = corr(anscombe(:,(i-1)*2+1),anscombe(:,(i-1)*2+2),'type','s');
    if i == 1
        labelt = 'accurate'
    elseif i == 2
        labelt = 'inflated'
    elseif i == 3
        labelt = 'deflated'
    elseif i == 4
        labelt = 'inflated'
    end
    title([ labelt ' r_p=' num2str(round(corr_p*1000)/1000) '; r_s=' num2str(round(corr_s*1000)/1000) ])
end
%% Figure 27.3
%NOP
%create random vector and showd the fisher transform

lots_of_corr_coefs = rand(1000,1)*2-1;
fisher_z_coefs = .5 * log( (1+lots_of_corr_coefs)./(1-lots_of_corr_coefs) );

figure
subplot(221)
hist(lots_of_corr_coefs,50)
xlabel('Correlation coefficient')
ylabel('Count')
axis square
set(gca,'xlim',[-1 1],'xtick',-1:.5:1)

subplot(222)
hist(fisher_z_coefs,50)
xlabel('Fisher-Z transformed coefficients')
ylabel('Count')
axis square
set(gca,'xlim',[-5 5],'xtick',-4:2:4)

subplot(223)
plot(lots_of_corr_coefs,fisher_z_coefs,'.')
xlabel('Correlation coefficient')
ylabel('Fisher-Z transformed coefficients')
set(gca,'xlim',[-1 1],'xtick',-1:.5:1)
axis square

subplot(224)
plot(atanh(lots_of_corr_coefs),fisher_z_coefs,'.')
xlabel('atanh')
ylabel('Fisher-Z')
r=corr(atanh(lots_of_corr_coefs),fisher_z_coefs);
legend([ 'Correlation = ' num2str(r) ])
axis square
set(gca,'xtick',-4:2:4,'ytick',-4:2:4)
axis([-4 4 -4 4])

end



function [convolution_result_fft1, convolution_result_fft2, r_pearson, r_spearman] = compute_convolution_onepair(sensor1, sensor2, EEG, half_of_wavelet_size, n_wavelet, n_data, n_convolution, wavelet_cycles, time,centerfreq, trial2plot,times2plot, draw_all_single_correlations)
%% Figure 27.4
%times2plot = dsearchn(EEG.times',[-300 1200]');
% global centerfreq;
% global trial2plot;
% global times2plot;
%times2plot =  dsearchn(EEG.times',[ EEG.times(1) EEG.times(end)]');
%trial2plot = EEG.trials;

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

function [ret_val] =  save_and_display_corr_mat(corrmat, initialrow, nbofrows, channel_labels, centerfreq, wavelet_cycles, subjandcond, typeofcorr)
%% Draws correlation matrix with contour using imagesc
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