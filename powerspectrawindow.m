%% Preprocesssing iEEG data from Natus acquisition system
% Author: Jaime Gomez Ramirez, Ph.D., mathematical neuroscience
% The Hospital for Sick Chldren, Toronto, ON
% email address: jd.gomezramirez@gmail.com  
% Website: aslab.org/~gomez
% November 2015; Last revision: 16-February-2015
function [] = powerspectrawindow()
%% powerspectrawindow Calculate the power spectra for selected (all) channels for a slideing window
% compare with powerspectra.m which is for the entire signal, no for
% sliding windows
% 1. Load initial data
% 2. Set N number of frequencies, we need only 100
% 3. Save mat files with the fourier coefficients and figures
%% Figure 11.5
%% 1. Load epoch and Initialize data
disp('Loading data....')
[myfullname, EEG_study, channel_labels] = initialize_EEG_variables();
EEG = EEG_study;
times2plot =  dsearchn(EEG.times',[ EEG.times(1) EEG.times(end)]');
trial2plot = EEG.trials;
%[pathstr, subjandcond, matext] = fileparts(myfullname);
disp('Data Loaded OK!');
%Choose which channels we want to calculate the power spectrum
initchan = 2;
chantouse = EEG.nbchan;
chantouse = 4;
% N= actual frequencies in Hz that will be returned by the
% Fourier transform. The number of unique frequencies we can measure is
% exactly 1/2 of the number of data points in the time series (plus DC).
N = EEG.pnts;
% sampling rate in Hz
srate   = EEG.srate; 
% Nyquist frequency -- the highest frequency you can measure in the data
nyquist = srate/2; 
frequencies = linspace(0,nyquist,N/2+1);
%vector of N points equally spaced when data is entire session, no windows
%time = ((1:N)-1)/N;
tot_secs = round(N/srate);

%% 2. Run short simul
% N = 60*srate; %10 seconds
N = 100;% number of frequencies to calculate
% We do not need Nyquist frequency , just from 0 to 50
snyq = 100;
%synqhalf = snyq/2+1;
synqhalf = snyq
%frequencies = linspace(0,snyq,snyq/2+1);
frequencies = linspace(0,snyq,snyq);

%% Loop for the specified channels, convolution for N frequencies
disp(['Calculating Fourier coefficients for channel=' num2str(initchan) ' to channel= ' num2str(chantouse)])
for chani=initchan:chantouse
    sensor2use = channel_labels(chani);     
    data_all = EEG.data(chani,:,trial2plot);
    fourier = zeros(size(data_all));   
    % Fourier transform is dot-product between sine wave and data at each frequency
    disp(['Calculating Fourier coefficients for channel:' sensor2use ' for Frequency=1, to Frequency=' num2str(N)])
    time = (1:srate)/srate;
    for fi=1:N
        disp(['Calcuating FT=sin*chani for freq=' num2str(fi) ' / ' num2str(N)]);
        sine_wave   = exp(-1i*2*pi*(fi-1).*time);
        % Build the window 
        fourierw = zeros(1,tot_secs);
        %size of the sliding windows in seconds
        winsize = 1;
        % time window selection, winsize = number of seconds
        for sec=1:winsize:tot_secs-1
               winini = (sec-1)*srate +1;
               winend = sec*srate;
               window = [winini winend];
               data = EEG.data(chani,winini:winend,trial2plot);
               partsum = sum(sine_wave.*data);
               fourierw(sec) = partsum;   
        end
        fourierw = sum(fourierw)/(tot_secs - 1);
        fourier(fi) = fourierw;
%       sine_wave   = exp(-1i*2*pi*(fi-1).*time);
%       fourier(fi) = sum(sine_wave.*data);
    end
    disp(['DONE with Fourier coefficients for channel:' sensor2use])
    
    %normalize fourier
    fourier=fourier/N;
    
    disp(['Saving power spetrum for channel=' num2str(chani) ' in ' myfullname])
    %Save object  powerspec as mat file
    [mat_pathstr,mat_name,mat_ext] = fileparts(myfullname);
    condsubjdate = strtok(mat_name, 'EEG_cut_');
    cd(mat_pathstr);
    if ~exist(mat_name, 'dir')
        mkdir(mat_name);
    end
    mat_pathstr =  fullfile(mat_pathstr,mat_name);
    cd(mat_pathstr);

    % abs fourier coeff
    it_angleft= angle(fourier(1:N/2+1));
    it_absft = abs(fourier(1:N/2+1)).^2;
    
    it_angleft= angle(fourier(1:synqhalf));
    it_absft = abs(fourier(1:synqhalf)).^2;
    
    mean_absft(chani-1) = mean(it_absft)
    std_absft(chani-1) = std(it_absft)
    mean_angleft(chani-1) = mean(it_angleft)
    std_angleft(chani-1) = std(it_angleft)
    disp(['Mean angle and abs =',  mean(it_angleft), std(it_angleft)])
    %% Save figures
    %displaypowerspec = 1 to display power spectrum charts, one per channel, if
    %0 only savemat file with data
    displaypowerspec = 1;
    if displaypowerspec == 1
        disp('Displaying the results')
        h = figure
        subplot(221)
        plot(data_all,'-o')
        set(gca,'xlim',[0 EEG.pnts])
        xlabel('Time')
        ylabel('Amplitude')
        title(['Time domain, Channel=' sensor2use])
        
        subplot(222)
        plot3(frequencies,it_angleft(1:synqhalf),it_absft(1:synqhalf),'-o','linew',3)
        grid on
        xlabel('Frequency (Hz)')
        ylabel('Phase')
        zlabel('Power')
        title('3-D of FT')
        view([20 20])
        
        subplot(223)
        bar(frequencies,it_absft(1:synqhalf))
        %set(gca,'xlim',[-5 105])
        set(gca,'xlim',[-5 70])
        xlabel('Frequency (Hz)')
        ylabel('Power')
        title('Power spectrum discrete FT')
        
        subplot(224)
        bar(frequencies,it_angleft(1:synqhalf))
        %set(gca,'xlim',[-5 105])
        set(gca,'xlim',[-5 70])
        xlabel('Frequency (Hz)')
        ylabel('Phase angle')
        set(gca,'ytick',-pi:pi/2:pi)
        title('Phase spectrum discrete FT')
        disp('Saving fourier coefficients and figure in figures')
        if ~exist('figures', 'dir')
            mkdir('figures');
        end
       
        olderfolder = cd('figures');
        matfh = sprintf('Powerspectrum_%s%s', sensor2use{1}, '.mat');
        disp(['Saving EEG epoch in file:  ', matfh])
        %save(matfh,'EEG', 'fourier') 
        matfigh = sprintf('Fig_Powerspectrum_%s%s', sensor2use{1}, '.fig')
        %savefig(h,matfigh);
        pause(5)
        close(h);
        cd(olderfolder); 
end
 
end
%disp(['Saving mat file:' powerspecfile ' in folder:' mat_pathstr]); 
%save(powerspecfile,'condsubjdate','absft','mean_absft','std_absft','angleft','mean_angleft','std_angleft')
%disp('DONE!')   
end