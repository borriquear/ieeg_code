%% Fourier transform analysis, to get power spectra and time frequency spectrum per
%  channel

function [] = fouriertransformanalysis()
%% fouriertransformanalysis : function that extracts power from the time series
% using the Fourier Transform

%% 1. Load epoch and Initialize data
disp('Loading the EEG data....')
[myfullname, EEG, channel_labels] = initialize_EEG_variables();
% condition, subject, date
% get condition, sesssion, patient and date
[eegpathname eegfilename eegextname]= fileparts(myfullname);
[eegcond, eegpatient,eegdate,eegsession ] = getsessionpatient(eegfilename)

%EEG = EEG_study;
times2plot =  dsearchn(EEG.times',[ EEG.times(1) EEG.times(end)]');
trial2plot = EEG.trials;
%[pathstr, subjandcond, matext] = fileparts(myfullname);
disp('Data Loaded OK!');
disp(['Total number of channels=' num2str(EEG.nbchan)]);
% number of trials
triali = 1;
% sampling rate
srate = EEG.srate;
disp(['Sampling rate=', num2str(srate)])
%total seconds of the epoch
tot_seconds = floor(EEG.pnts/srate);
disp(['Total seconds of the session=' num2str(tot_seconds)])
% remaining points beyond tot_seconds
remainingpoints = EEG.pnts - tot_seconds*srate;
%time vector normalized
t = 0:1/srate:(tot_seconds + ((remainingpoints/srate)) );
% end -1 to have samne length as signal
t = t(1:end-1);
% t =horzcat(t,restsecondslist);
n = length(t); % = EEG.pnts
nyquistfreq = srate/2;
% There is one half only of frequencies of fourier coefficients
hz = linspace(0,nyquistfreq,floor(n/2)+1);

% Choose the channel(s) we want to analyze
chani = 2;
chanend = EEG.nbchan;
chanend = chani +2;
%chanend = chani;
tot_channels = chanend - chani+1;
%vector of channels (0,1) 1 when to display that channel
vectorofchannelsprint = zeros(1,EEG.nbchan);
vectorofchannelsprint(chani:chanend) = 1;
tot_rows = tot_channels;
% figure all channel bands
allchbands = figure;
hbars = figure;
for irow =chani:chanend
    hchann(chani) = figure;
    signal = EEG.data(chani,:,triali);
    signalX = zeros(size(signal));
    % fourier time is not the time in seconds but the time normalized
    fouriertime = (0:n)/n;
    fouriertime = fouriertime(1:end-1);
    % discrete fourier Transform
    disp(['Calculating the DFT for channel=' num2str(chani)])
    disp('Reminder: For long time series use FFT instead')
    %fftidx =1 use fft, much faster
    fftidx  = 1;
    if fftidx == 1
        %% FFT fft(X) computes the discrete Fourier transform (DFT) of X using a fast Fourier transform (FFT) algorithm.
        if vectorofchannelsprint(chani) == 1
            disp('Calculating FFT...')
            signalXF = fft(signal)/n;
            figure(hchann(chani))
            title([' Condition=' eegcond ', Patient =' eegpatient ', Date' eegdate ', Session' eegsession ] );
            %subplot(tot_rows,4,1 + (4*(irow-1)))
            % plot the signal time series
            subplot(2,1,1)
            plot(t,signal)
            xlabel('Time (s)'), ylabel('Amplitude')
            %set(gca,'ylim',[-0.2 6.0])
            set(gca,'xlim',[0 tot_seconds])
            legend({'time series'})
            %plot fourier coefficients
            %subplot(tot_rows,4,2 + (4*(irow-1)))
            subplot(2,1,2)
            plot(hz,2*abs(signalXF(1:length(hz))),'r')
            xlabel('Frequencies (Hz)'), ylabel('Amplitude')
            set(gca,'xlim',[0 50])
            legend({'fast Fourier transform'})
            msgtitle = sprintf('Correlation Coefficient between pairs: %s ', typeofcorr);
            title(msgtitle);
        end
    else
        % not fast fourier very inefficient just to see how conceptually
        % works
        for fi = 1:length(signalX)
            csw=exp(-1i*2*pi*(fi-1)*fouriertime);
            signalX(fi) = sum(csw.*signal)/n;
            formatSpec = 'Calculating freq=%d/totalfreqs=%d';
            sprintf(formatSpec, fi, EEG.pnts)
        end
    end
    %% Multitaper method for power spectra
    % multitaper method may be useful in situations of low SNR
    %The multitaper method is an extension of the Fourier transform, in which the Fourier transform is computed several times, each time tapering the data using a different taper. The tapers are taken from the "Slepian sequence," in which each taper is orthogonal to each other taper.
    % how much noise to add
    %it estimates spectral peaks that are wider than the frequencies present in the original time series (this is called spectral leakage or frequency smearing). Thus, if isolating closely spaced frequencies is important, the multitaper method may not be a preferable option.
    %noisefactor = 20;
    % calculate power spectra with multitaper method
    multitaper = 0;
    if multitaper == 1
        disp(['Multitaper method to calculate the power spectrum for channel=' num2str(irow)])
        hmt(chani) = figure;
        % define Slepian tapers.
        tapers = dpss(n,3)';
        % initialize multitaper power matrix
        mtPow = zeros(floor(n/2)+1,1);
        hz = linspace(0,srate/2,floor(n/2)+1);
        % loop through tapers
        for tapi = 1:size(tapers,1)-1 % -1 because the last taper is typically not used
            % scale the taper for interpretable FFT result
            temptaper = tapers(tapi,:)./max(tapers(tapi,:));
            
            % FFT of tapered data
            x = abs(fft(signal.*temptaper)/n).^2;
            
            % add this spectral estimate to the total power estimate
            mtPow = mtPow + x(1:length(hz))';
        end
        % Because the power spectra were summed over many tapers,
        % divide by the number of tapers to get the average.
        mtPow = mtPow./tapi;
        
        % now compute the 'normal' power spectra using one taper
        hann   = .5*(1-cos(2*pi*(1:n)/(n-1)));
        x      = abs(fft(signal.*hann)/n).^2;
        regPow = x(1:length(hz)); % regPow = regular power
        
        % Now plot both power spectra. Note that power is plotted here instead of
        % amplitude because of the amount of noise. Try plotting amplitude by
        % multiplying the FFT result by 2 instead of squaring.
        clf
        plot(hz,mtPow,'.-'), hold on
        plot(hz,regPow,'r.-')
        set(gca,'xlim',[0 50])
        legend({'single taper';'multi taper'})
    end
    
    %% Extract information about specific frequencies.
    % frequencies vector
    %figure
    f = 0:50;
    frex_idx = sort(dsearchn(hz',f'));
    requested_frequences = 2*abs(signalXF(frex_idx));
    
    %clf
    figure(hbars);
    % subplot(tot_rows,3,3*irow)
    % %subplot(111)
    bar(requested_frequences)
    xlabel('Frequencies (Hz)'), ylabel('Amplitude')
    freq_bands = ['t' 'd' 'a' 'b' 'g' ]
    set(gca,'xtick',1:length(frex_idx),'xticklabel',cellstr(num2str(round(hz(frex_idx))')))
    % %hold on
    x1 = 4;
    x2 = 8;
    x3 = 12;
    x4 = 40;
    y1=get(gca,'ylim');
    hold on
    plot([x1 x1],y1)
    plot([x2 x2],y1)
    plot([x3 x3],y1)
    plot([x4 x4],y1)
    
    %% Time - frequency analysis, The short-time Fourier transform
    sft = 0;
    if sft == 1
        % define 'chunks' of time for different frequency components
        timechunks = round(linspace(1,length(t),length(f)+1));
        %figure
        %subplot(211)
        %subplot(tot_rows,4,3*irow)
        %plot(t,signal)
        %xlabel('Time (s)'), ylabel('Amplitude')
        
        % short-time FFT parameters
        fftWidth_ms = 1000;
        fftWidth    = round(fftWidth_ms/(1000/srate)/2);
        Ntimesteps  = tot_seconds; % number of time widths
        centertimes = round(linspace(fftWidth+1,length(t)-fftWidth,Ntimesteps));
        
        % frequencies in Hz
        hz = linspace(0,srate/2,fftWidth-1);
        
        % initialize matrix to store time-frequency results
        tf = zeros(length(hz),length(centertimes));
        
        % Hann window to taper the data in each window for the FFT
        hanwin = .5*(1-cos(2*pi*(1:fftWidth*2)/(fftWidth*2-1)));
        
        % loop through center time points and compute Fourier transform
        for ti=1:length(centertimes)
            % get data from this time window
            %temp = data(centertimes(ti)-fftWidth:centertimes(ti)+fftWidth-1);
            temp = EEG.data(chani,centertimes(ti)-fftWidth:centertimes(ti)+fftWidth-1,triali);
            % Fourier transform
            x = fft(hanwin.*temp)/fftWidth*2;
            
            % enter amplitude into tf matrix
            tf(:,ti) = 2*abs(x(1:length(hz)));
        end
        % plot the time-frequency result
        %subplot(212)
        %subplot(tot_rows,4,3*irow)
        hsft(chani) = figure;
        
        contourf(t(centertimes),hz,tf,1)
        set(gca,'ylim',[0 max(t(centertimes))],'clim',[0 1],'xlim',[0 tot_seconds])
        xlabel('Time (s)'), ylabel('Frequency (Hz)')
        legend('sFT')
        
        % The figures in the book use a reverse colorscaling,
        % which can be obtained using the following code.
        c = gray; % 'gray' is the colormap
        colormap(c(end:-1:1,:)) % reversed colormap
    end
    %end short fourier
    %% Morlet Wavelet
    %Chapter 6.2, Figure 6.5
    % create a linear chirp
    %srate   = 1000;
    %t       = 0:1/srate:6;
    %f   = [2 8]; % frequencies in Hz
    %chirpTS = sin(2*pi.*linspace(f(1),mean(f),length(t)).*t);
    
    % The wavelet has its own time vector, but
    % the sampling rate MUST be the same as for the signal.
    % for a wave between -2 and 2s
    timewvi = -2;
    timewve = 2;
    wavetime = timewvi:1/srate:timewve;
    
    % width of the Gaussian that tapers the sine wave
    wvcycles = 4;
    wvfreq = 8;
    w = 2*( wvcycles/(2*pi*wvfreq) )^2;
    
    % create the complex Morlet wavelet
    cmw = exp(1i*2*pi*wvfreq.*wavetime) .* exp( (-wavetime.^2)/w );
    
    % half of the length of the wavelet
    halfwavsize = floor(length(wavetime)/2);
    
    % Perform time-domain convolution...
    % the signal must be zero-padded
    signalpad = [zeros(1,halfwavsize) signal zeros(1,halfwavsize)];
    %chirpTSpad = [zeros(1,halfwavsize) chirpTS zeros(1,halfwavsize)];
    % initialize results matrix
    %convres = zeros(size(chirpTSpad));
    convres = zeros(size(signalpad));
    
    % and now the actual convolution
    for i=halfwavsize+1:length(signal)+halfwavsize-1
        convres(i) = sum( signalpad(i-halfwavsize:i+halfwavsize) .* cmw );
    end
    % the final step of convolution is to trim the edges
    convres = convres(halfwavsize:end-halfwavsize-1);
    
    % and plot the results
    % %clf
    % subplot(211)
    % plot(t,signal)
    % % the amplitude of the result of the convolution of the signal and a
    % % wvfreq=5 wavelet
    % subplot(212)
    % plot(t,abs(convres))
    % xlabel('Time (s)'), ylabel('Amplitude (unscaled)')
    % for th vonvolution theorem we can multiply the frequency domain and is
    % equibalent to the dot product in the time domain
    %% To normalize the amplitude
    
    Lconv = length(t)+length(wavetime)-1;
    
    % Fourier spectra of the two signals
    cmwX = fft(cmw,Lconv);
    cmwX = cmwX./max(cmwX);
    
    % and their multiplied inverse
    convres4 = ifft( fft(signal,Lconv).*cmwX );
    convres4 = convres4(halfwavsize:end-halfwavsize-1);
    
    %clf
    %subplot(tot_rows,4,1 + (4*(irow-1))), hold on
    %figure
    
    %plot(t,signal), hold on
    if vectorofchannelsprint(chani) == 1
        figure(hchann(chani));
        subplot(2,1,1), hold on
        plot(t,2*abs(convres4),'r')
    end
    %% time frequency analysis with Morlet wavelet
    disp('Time frequency analysis with Morlet wavelet...')
    % number and range of frequencies for analysis
    nfrex = 50;
    frex  = logspace(log10(2),log10(nfrex),nfrex);
    
    % initialize...
    tf = zeros(nfrex,length(signal));
    
    % Fourier spectrum of the signal. Note that because this does not change as
    % a function of wavelet frequency, it needs to be computed only once.
    signalx = fft(signal,Lconv);
    
    % loop through frequencies
    for fi=1:nfrex
        % compute normalized Fourier spectra of wavelet
        disp(['compute normalized Fourier spectra of wavelet fi=' num2str(fi) ' / ' num2str(nfrex)])
        w = 2*( 5/(2*pi*frex(fi)) )^2;
        cmwX = fft(exp(1i*2*pi*frex(fi).*wavetime) .* exp( (-wavetime.^2)/w ), Lconv);
        cmwX = cmwX./max(cmwX);
        % convolution and extract amplitude
        convres = ifft( signalx .* cmwX );
        tf(fi,:) = 2*abs(convres(halfwavsize:end-halfwavsize-1));
    end
    disp('end loop to calculate tf')
    % w=600; h=400;
    % set(figure(4),'Position',[400 60 w h])
    %subplot(tot_rows,4,4  + (4*(irow-1)))
    htimefre(chani)= figure;
    disp('calculating contour...')
    contourf(t,frex,tf,40,'linecolor','none')
    set(gca,'clim',[0 1]), colorbar
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    figure(allchbands);
    hold on
    plot(hz,2*abs(signalXF(1:length(hz))),'r')
    xlabel('Frequencies (Hz)'), ylabel('Amplitude')
    set(gca,'xlim',[0 50])
    legend({'fast Fourier transform'})
    disp(['DONE with channel=' num2str(irow)])
    totalsignalXF = sum(2*abs(signalXF(1:length(hz))));
end
%end of channels
figure(allchbands);
hamperband = 0;
x1 = 4;
x2 = 8;
x3 = 12;
x4 = 40;
y1=get(gca,'ylim')
hold on
plot([x1 x1],y1)
plot([x2 x2],y1)
plot([x3 x3],y1)
plot([x4 x4],y1)
% deltaband = mean(2*abs(totalsignalXF(x1:x2)))
% alphaband = mean(2*abs(totalsignalXF(x2:x3)))
% betaband = mean(2*abs(totalsignalXF(x3:x4)))
% gammaband= mean(2*abs(totalsignalXF(x4:end)))
end

function [eegcond, eegpatient,eegdate,eegsession ] = getsessionpatient(eegfilename)
%%     [eegcond, eegpatient,eegdate,eegsession ] = getsessionpatient(eegfilename)  get patient name etc from file name

[eegbl remain] = strtok(eegfilename, 'EEG_cut_');
[eegblcond remain]= strtok(remain, '_');
[eegblcond2 remain]= strtok(remain, '_');
if (strcmpi(eegblcond2,'PRE') == 1) || (strcmpi(eegblcond2,'POST') == 1)
    eegcond = strcat(eegbl, eegblcond);
    eegcond = strcat(eegcond, eegblcond2);
    [eegpatient remain]= strtok(remain, '_');
    [eegdate remain]= strtok(remain, '_');
    eegsession = strtok(remain, '_');
else
    eegcond = strcat(eegbl, eegblcond);
    eegpatient = eegblcond2;
    [eegdate eegsession]= strtok(remain, '_');
    eegsession = strtok(eegsession, '_');
end
end