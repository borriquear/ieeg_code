%Scripts for powerbased connectivity with iEEG data
function [ ] =  powerbasedconnectivity()

%% Figure 27.2
%Distribution of power and logPower per channel overtime 
% Display only channels that have a power or log power normal distribution

%load data
%mydir = 'C:\Users\shelagh\Desktop\patients\Cordeiro, Jessica\Session 1 day 1 20 October\dp_cj28_20151020_s1'
[myfullname, EEG_study, channel_labels, initchan, centerfreq, chantouse, time] = initialize_EEG_variables();
[half_of_wavelet_size, n_wavelet, n_data n_convolution, wavelet_cycles] = initialize_wavelet(time, EEG_study);

for chani=initchan:chantouse   
    %sensor2use = 'LMT3'
    sensor2use = channel_labels(chani)
    % setup wavelet convolution and outputs   
    % FFT of data (note: this doesn't change on frequency iteration)
    fft_data = fft(reshape(EEG_study.data(strcmpi(sensor2use,{EEG_study.chanlocs.labels}),:,:),1,EEG_study.pnts*EEG_study.trials),n_convolution);   
    % create wavelet and run convolution
    fft_wavelet            = fft(exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq))^2)),n_convolution);
    convolution_result_fft = ifft(fft_wavelet.*fft_data,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    convolution_result_fft = abs(reshape(convolution_result_fft,EEG_study.pnts,EEG_study.trials)).^2;
    
    % trim edges so the distribution is not driven by edge artifact outliers
    % (note: here we just use visual inspection to remove edges)
    convolution_result_fft = convolution_result_fft(100:end-100,:);
    
    % test for normal distribution, if you have the stats toolbox
    if exist('kstest','file')
        [h,p1] = kstest(convolution_result_fft(:));
        [h,p2] = kstest(log10(convolution_result_fft(:)));
        [h,p3] = kstest(randn(numel(convolution_result_fft),1));
        disp([ 'KS test for normality of power: '        num2str(p1) ' (>.05 means normal distribution) ' ])
        disp([ 'KS test for normality of log10(power): ' num2str(p2) ' (>.05 means normal distribution) ' ])
        disp([ 'KS test for normality of random data: '  num2str(p3) ' (>.05 means normal distribution) ' ])
    end
    if (p1 + p2 > 0.05)
        % plot distirbution of power data
        figure
        %subplot(chantouse-1,2,(chani-1)*2-1)
        subplot(1,2,1)
        hist(convolution_result_fft(:),100)
        title('Distribution of power values')
        axis square
        %set(gca,'xlim',[EEG.times(1) EEG.times(end)],'ylim',[min(real_data(:)) max(real_data(:))])
        xlabel('Power (\muV))')
        ylabel('Count')
        %subplot(chantouse-1,2,(chani-1)*2)
        subplot(1,2,2)
        hist(log10(convolution_result_fft(:)),100)
        title('Distribution of log_1_0power values')
        xlabel('log 10 Power')
        ylabel('Count')
        axis square
end
end 

%% Figure 27.4
%calculate channel correlations

%centerfreq = 10; % in Hz
trial2plot = 1;
times2plot = dsearchn(EEG_study.times',[EEG_study.times(1) EEG_study.times(end)]');
chantouse = 3;
initchan = 2;
convolution_result_fft = [];
convolution_result_fft1 = [];
convolution_result_fft2 = [];
for chani=initchan:chantouse
    for chanj=chani+1:chantouse
        sensor1 = channel_labels(chani);
        sensor2 = channel_labels(chanj);
       
        %esto no hace falta
        %times2plot = dsearchn(EEG_study.times',[-300 1200]');
        times2plot = dsearchn(EEG_study.times',[EEG_study.times(1) EEG_study.times(end)]');
        
        fft_data1 = fft(reshape(EEG_study.data(strcmpi(sensor1,{EEG_study.chanlocs.labels}),:,:),1,EEG_study.pnts*EEG_study.trials),n_convolution);
        fft_data2 = fft(reshape(EEG_study.data(strcmpi(sensor2,{EEG_study.chanlocs.labels}),:,:),1,EEG_study.pnts*EEG_study.trials),n_convolution);
        
        % create wavelet and run convolution
        fft_wavelet             = fft(exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq))^2)),n_convolution);
        convolution_result_fft  = ifft(fft_wavelet.*fft_data1,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq));
        convolution_result_fft  = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
        convolution_result_fft1 = reshape(convolution_result_fft,EEG_study.pnts,EEG_study.trials);
        
        fft_wavelet             = fft(exp(2*1i*pi*centerfreq.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq))^2)),n_convolution);
        convolution_result_fft  = ifft(fft_wavelet.*fft_data2,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq));
        convolution_result_fft  = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
        convolution_result_fft2 = reshape(convolution_result_fft,EEG_study.pnts,EEG_study.trials);
        
        % keep only requested time regions
        convolution_result_fft1 = convolution_result_fft1(times2plot(1):times2plot(2),:);
        convolution_result_fft2 = convolution_result_fft2(times2plot(1):times2plot(2),:);
        
        figure
        subplot(211)
        plot(EEG_study.times(times2plot(1):times2plot(2)),abs(convolution_result_fft1(:,trial2plot)).^2)
        hold on
        plot(EEG_study.times(times2plot(1):times2plot(2)),abs(convolution_result_fft2(:,trial2plot)).^2,'r')
        xlabel('Time (ms)')
        ylabel('Power')
        set(gca,'xlim',EEG_study.times(times2plot))
        legend([sensor1 sensor2] )
        title(['Power time series at freq=' num2str(centerfreq) ])
        
        subplot(223)
        plot(abs(convolution_result_fft1(:,trial2plot)).^2,abs(convolution_result_fft2(:,trial2plot)).^2,'.')
        xlabel([ sensor1 ])
        ylabel([ sensor2 ])
        r=corr(abs(convolution_result_fft1(:,trial2plot)).^2,abs(convolution_result_fft2(:,trial2plot)).^2,'type','p');
        %set(gca,'xlim',EEG_study.times(times2plot))
        title(['Pearson '  num2str(r) ])
        legend([ 'Pearson R = ' num2str(r) ]);
       
        subplot(224)
        plot(tiedrank(abs(convolution_result_fft1(:,trial2plot)).^2),tiedrank(abs(convolution_result_fft2(:,trial2plot)).^2),'.')
        xlabel([ sensor1 ])
        ylabel([ sensor2 ])
        r=corr(abs(convolution_result_fft1(:,trial2plot)).^2,abs(convolution_result_fft2(:,trial2plot)).^2,'type','s');
        legend([ 'Spearman Rho = ' num2str(r) ]);
        title(['Spearman '  num2str(r) ])
        %set(gca,'ylim',get(gca,'xlim'))
        
    end
end

%% Figure 27.5

% Compute how many time points are in one cycle, and limit xcov to this lag
nlags = round(EEG_study.srate/centerfreq);

% note that data are first tiedrank'ed, which results in a Spearman rho
% instead of a Pearson r. 
[corrvals,corrlags] = xcov(tiedrank(abs(convolution_result_fft1(:,trial2plot)).^2),tiedrank(abs(convolution_result_fft2(:,trial2plot)).^2),nlags,'coeff');

% convert correlation lags from indices to time in ms
corrlags = corrlags * 1000/EEG_study.srate; 

figure
plot(corrlags,corrvals,'-o','markerface','w')
hold on
plot([0 0],get(gca,'ylim'))
xlabel([   sensor1{1} ' leads --- Time lag in ms --- ' sensor2{1} ])
ylabel('Correlation coefficient')
title('Corr. coef. as a function of time lag')

%% Figure 27.6

%sensor1 = 'poz';
%sensor2 = 'fz';

% timewin1 = [ -300 -100 ]; % in ms relative to stim onset
% timewin2 = [  200  400 ];
% 
% centerfreq1 =  6; % in Hz
% centerfreq2 =  6;
% 
% 
% % convert time from ms to index
% timeidx1 = zeros(size(timewin1));
% timeidx2 = zeros(size(timewin2));
% for i=1:2
%     [junk,timeidx1(i)] = min(abs(EEG.times-timewin1(i)));
%     [junk,timeidx2(i)] = min(abs(EEG.times-timewin2(i)));
% end
% 
% % setup wavelet convolution and outputs
% time = -1:1/EEG.srate:1;
% half_of_wavelet_size = (length(time)-1)/2;
% 
% % FFT parameters
% n_wavelet     = length(time);
% n_data        = EEG.pnts*EEG.trials;
% n_convolution = n_wavelet+n_data-1;
% wavelet_cycles= 4.5;
% 
% % FFT of data (note: this doesn't change on frequency iteration)
% fft_data1 = fft(reshape(EEG.data(strcmpi(sensor1,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);
% fft_data2 = fft(reshape(EEG.data(strcmpi(sensor2,{EEG.chanlocs.labels}),:,:),1,EEG.pnts*EEG.trials),n_convolution);
% 
% % initialize output time-frequency data
% corrdata = zeros(EEG.trials,2);
% 
% % create wavelet and run convolution
% fft_wavelet            = fft(exp(2*1i*pi*centerfreq1.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq1))^2)),n_convolution);
% convolution_result_fft = ifft(fft_wavelet.*fft_data1,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq1));
% convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
% convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
% analyticsignal1        = abs(convolution_result_fft).^2;
% 
% fft_wavelet            = fft(exp(2*1i*pi*centerfreq2.*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*centerfreq2))^2)),n_convolution);
% convolution_result_fft = ifft(fft_wavelet.*fft_data2,n_convolution) * sqrt(wavelet_cycles /(2*pi*centerfreq2));
% convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
% convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
% analyticsignal2        = abs(convolution_result_fft).^2;
% 
% 
% % Panel A: correlation in a specified window
% tfwindowdata1 = mean(analyticsignal1(timeidx1(1):timeidx1(2),:),1);
% tfwindowdata2 = mean(analyticsignal2(timeidx2(1):timeidx2(2),:),1);
% figure
% subplot(121)
% plot(tfwindowdata1,tfwindowdata2,'.')
% axis square
% title([ 'TF window correlation, r_p=' num2str(corr(tfwindowdata1',tfwindowdata2','type','p')) ])
% xlabel([ sensor1 ': ' num2str(timewin1(1)) '-' num2str(timewin1(2)) '; ' num2str(centerfreq1) ' Hz' ])
% ylabel([ sensor2 ': ' num2str(timewin2(1)) '-' num2str(timewin2(2)) '; ' num2str(centerfreq2) ' Hz' ])
% 
% 
% % also plot rank-transformed data
% subplot(122)
% plot(tiedrank(tfwindowdata1),tiedrank(tfwindowdata2),'.')
% axis square
% xlabel([ sensor1 ': ' num2str(timewin1(1)) '-' num2str(timewin1(2)) '; ' num2str(centerfreq1) ' Hz' ])
% ylabel([ sensor2 ': ' num2str(timewin2(1)) '-' num2str(timewin2(2)) '; ' num2str(centerfreq2) ' Hz' ])
% title([ 'TF window correlation, r_p=' num2str(corr(tfwindowdata1',tfwindowdata2','type','s')) ])
% 
% 
% % panel B: correlation over time
% corr_ts = zeros(size(EEG.times));
% for ti=1:EEG.pnts
%     corr_ts(ti) = corr(analyticsignal1(ti,:)',analyticsignal2(ti,:)','type','s');
% end
% 
% figure
% plot(EEG.times,corr_ts)
% set(gca,'xlim',[-200 1200])
% xlabel('Time (ms)'), ylabel('Spearman''s rho')
% 
% 
% 
% % Panel C: exploratory time-frequency power correlations
% times2save = -200:25:1200;
% frex = logspace(log10(2),log10(40),20);
% 
% 
% times2save_idx = zeros(size(times2save));
% for i=1:length(times2save)
%     [junk,times2save_idx(i)] = min(abs(EEG.times-times2save(i)));
% end
% 
% % rank-transforming the data can happen outside the frequency loop
% seeddata_rank = tiedrank(tfwindowdata2);
% 
% % initialize output correlation matrix
% expl_corrs = zeros(length(frex),length(times2save));
% 
% for fi=1:length(frex)
%     
%     % get power (via wavelet convolution) from signal1
%     fft_wavelet            = fft(exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frex(fi)))^2)),n_convolution);
%     convolution_result_fft = ifft(fft_wavelet.*fft_data1,n_convolution) * sqrt(wavelet_cycles /(2*pi*frex(fi)));
%     convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
%     convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
%     analyticsignal1        = abs(convolution_result_fft).^2;
%     
%     for ti=1:length(times2save)
%         expl_corrs(fi,ti) = 1-6*sum((seeddata_rank-tiedrank(analyticsignal1(times2save_idx(ti),:))).^2)/(EEG.trials*(EEG.trials^2-1));
%     end
% end
% 
% figure
% contourf(times2save,frex,expl_corrs,40,'linecolor','none')
% set(gca,'clim',[-.4 .4],'yscale','log','ytick',round(logspace(log10(frex(1)),log10(frex(end)),8)))
% xlabel('Time (ms)'), ylabel('Frequency (Hz)')
% title([ 'Correlation over trials from seed ' sensor2 ', ' num2str(centerfreq2) ' Hz and ' num2str(timewin2(1)) '-' num2str(timewin2(2)) ' ms' ])
% colorbar

%% Figure 27.7
% Time-frequency map of partial correlation coefficients between power at two electrodes, seed_chan and target chan
%(panel A), while holding
% (panel B), corr. coeff holding a third electrode constant (control_chan) 
% Partial correlations can be used to test network-level hypotheses involving
% more than two electrodes, or they can be used to minimize the effects of volume conduction on powerbased
% connectivity by removing the variance between two electrodes (in this case,
% seed-target chans) that is shared with — and likely volume-conducted from — a neighboring electrode (control_chan).
seed_chan     = 'LHD1';
target_chan   = 'LHD2';
control_chan  = 'LHD3';

clim = [0 .6];

% wavelet parameters
min_freq = 2;
max_freq = 40;
num_frex = 15;

% downsampled times
%times2save = -200:50:800;
%times2save = EEG_study.times; % uncomment this line for figure 27.8
totalsecsepoch = round(times2plot(2)/EEG_study.srate);
init_time = times2plot(1); %-1; %times2plot(1)
end_time = (totalsecsepoch/10)*EEG_study.srate; %only 10% of the entire epoch %times2plot(2)
times2save = init_time:EEG_study.srate/10:end_time;
times2save = EEG_study.times;
% other wavelet parameters
frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);

%time = init_time:1/EEG_study.srate:end_time;
time = EEG_study.times;
time= 0:1:1000;
%time = times2save;
%half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters
% n_wavelet     = length(time);
% n_data        = EEG.pnts*EEG.trials;
% n_convolution = n_wavelet+n_data-1;
% wavelet_cycles= 4.5;

[half_of_wavelet_size, n_wavelet, n_data n_convolution, wavelet_cycles] = initialize_wavelet(time, EEG_study);

times2saveidx = dsearchn(EEG_study.times',times2save');

% FFT of data (note: this doesn't change on frequency iteration)
fft_data_seed = fft(reshape(EEG_study.data(strcmpi(seed_chan,{EEG_study.chanlocs.labels}),:,:),1,EEG_study.pnts*EEG_study.trials),n_convolution);
fft_data_trgt = fft(reshape(EEG_study.data(strcmpi(target_chan,{EEG_study.chanlocs.labels}),:,:),1,EEG_study.pnts*EEG_study.trials),n_convolution);
fft_data_ctrl = fft(reshape(EEG_study.data(strcmpi(control_chan,{EEG_study.chanlocs.labels}),:,:),1,EEG_study.pnts*EEG_study.trials),n_convolution);


% initialize output time-frequency data
tf_corrdata = zeros(length(frequencies),length(times2save),2);

for fi=1:length(frequencies)
    
    % create wavelet and get its FFT
    fft_wavelet = fft(exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi),n_convolution);
    
    % convolution for all three sites (save only power)
    convolution_result_fft = ifft(fft_wavelet.*fft_data_seed,n_convolution) * sqrt(wavelet_cycles /(2*pi*frequencies(fi)));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    conv_result_seed = abs(reshape(convolution_result_fft,EEG_study.pnts,EEG_study.trials)).^2;
    
    convolution_result_fft = ifft(fft_wavelet.*fft_data_trgt,n_convolution) * sqrt(wavelet_cycles /(2*pi*frequencies(fi)));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    conv_result_trgt = abs(reshape(convolution_result_fft,EEG_study.pnts,EEG_study.trials)).^2;
    
    convolution_result_fft = ifft(fft_wavelet.*fft_data_ctrl,n_convolution) * sqrt(wavelet_cycles /(2*pi*frequencies(fi)));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
    conv_result_ctrl = abs(reshape(convolution_result_fft,EEG_study.pnts,EEG_study.trials)).^2;
    
    % downsample and rank transform all data
    conv_result_seed = tiedrank(conv_result_seed(times2saveidx,:)')';
    conv_result_trgt = tiedrank(conv_result_trgt(times2saveidx,:)')';
    conv_result_ctrl = tiedrank(conv_result_ctrl(times2saveidx,:)')';
    
    for ti=1:length(times2save)
        
        % compute bivariate correlations
        r_st = 1-6*sum((conv_result_seed(ti,:)-conv_result_trgt(ti,:)).^2);
        r_sc = 1-6*sum((conv_result_seed(ti,:)-conv_result_ctrl(ti,:)).^2);
        r_tc = 1-6*sum((conv_result_ctrl(ti,:)-conv_result_trgt(ti,:)).^2);
        
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
    set(gca,'clim',clim,'xlim',[times2save(1) times2save(end)],'yscale','log','ytick',logspace(log10(frequencies(1)),log10(frequencies(end)),6),'yticklabel',round(logspace(log10(frequencies(1)),log10(frequencies(end)),6)*10)/10)
    %set(gca,'clim',clim,'xlim',[init_time end_time],'yscale','log','ytick',logspace(log10(frequencies(1)),log10(frequencies(end)),6),'yticklabel',round(logspace(log10(frequencies(1)),log10(frequencies(end)),6)*10)/10)

    axis square
    if i==1
        title([ 'Correlation ' seed_chan ' , ' target_chan ])
    else
        title([ 'Conditional Partial correlation ' seed_chan ', ' target_chan '| ' control_chan ])
    end
    xlabel('Time (ms)'), ylabel('Frequency (Hz)')
end

%% Figure 27.8

% Re-run the code for the previous figure but comment out the 
% following line towards the top:
% times2save = EEG.times; % uncomment this line for figure 27.8
% Then run this section of code.


ds_timesidx = dsearchn(EEG_study.times',(-200:50:800)'); % ds = downsampled
ds_timesidx = dsearchn(EEG_study.times',(0:50:8000)');
[~,lofreq]  = min(abs(frequencies-4.7));
[~,hifreq]  = min(abs(frequencies-32));

figure

subplot(221)
contourf(times2save,frequencies,squeeze(tf_corrdata(:,:,2)),40,'linecolor','none')
hold on
plot(get(gca,'xlim'),frequencies([lofreq lofreq]),'k--')
plot(get(gca,'xlim'),frequencies([hifreq hifreq]),'k--')
set(gca,'clim',clim,'xlim',[0 8000],'yscale','log','ytick',logspace(log10(frequencies(1)),log10(frequencies(end)),6),'yticklabel',round(logspace(log10(frequencies(1)),log10(frequencies(end)),6)*10)/10)
title('Original (256 Hz)')


subplot(222)
contourf(times2save(ds_timesidx),frequencies,squeeze(tf_corrdata(:,ds_timesidx,2)),40,'linecolor','none')
hold on
plot(get(gca,'xlim'),frequencies([lofreq lofreq]),'k--')
plot(get(gca,'xlim'),frequencies([hifreq hifreq]),'k--')
set(gca,'clim',clim,'xlim',[-0 8000],'yscale','log','ytick',logspace(log10(frequencies(1)),log10(frequencies(end)),6),'yticklabel',round(logspace(log10(frequencies(1)),log10(frequencies(end)),6)*10)/10)
title('Down-sampled (20 Hz)')


subplot(223)
plot(EEG_study.times,squeeze(tf_corrdata(lofreq,:,2)))
hold on
plot(EEG_study.times(ds_timesidx),squeeze(tf_corrdata(lofreq,ds_timesidx,2)),'ro-','markerface','w')
title('Effect of downsampling on low-frequency activity')
set(gca,'xlim',[0 8000],'ylim',[.25 .65])


subplot(224)
plot(EEG_study.times,squeeze(tf_corrdata(hifreq,:,2)))
hold on
plot(EEG_study.times(ds_timesidx),squeeze(tf_corrdata(hifreq,ds_timesidx,2)),'ro-','markerface','w')
title('Effect of downsampling on high-frequency activity')
set(gca,'xlim',[0 8000],'ylim',[-.1 .6])
legend({'Original (1000 Hz)';'Down-sampled (20 Hz)'})

%% end.
