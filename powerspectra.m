function [] = powerspectra()
%% Calculate the power spectra for selected (all) channels
%% Figure 11.5
%% 1. Load epoch and Initialize data
disp('Loading data....')
[myfullname, EEG_study, channel_labels] = initialize_EEG_variables();
EEG = EEG_study;
times2plot =  dsearchn(EEG.times',[ EEG.times(1) EEG.times(end)]');
trial2plot = EEG.trials;
%[pathstr, subjandcond, matext] = fileparts(myfullname);
disp('Data Loaded OK!');

initchan = 6;
chantouse = EEG.nbchan;
chantouse = initchan;
N = EEG.pnts;% length of sequence    
time = ((1:N)-1)/N;
%time = ((1:N)-1);
srate   = EEG.srate;        % sampling rate in Hz
nyquist = srate/2;    % Nyquist frequency -- the highest frequency you can measure in the data
%partition = 10;
%nyquist_short = srate/partition;
frequencies = linspace(0,nyquist,N/2+1);
%frequencies_short = linspace(0,nyquist_short,N/partition+1);
endidx = round(N/2) + 1;
% power vector  fourier coeff, mean and std
absft = zeros(EEG.nbchan-1,endidx);
mean_absft = zeros(1, EEG.nbchan-1);
std_absft = zeros(1, EEG.nbchan-1);
%angle fourier coeff
angleft = zeros(EEG.nbchan-1,endidx);
mean_angleft = zeros(1, EEG.nbchan-1);
std_angleft = zeros(1, EEG.nbchan-1);
disp(['Calculating Fourier coefficients for channel=' num2str(initchan) ' to channel= ' num2str(chantouse)])
for chani=initchan:chantouse
%disp(['.....Channel=' sensor2use])
sensor2use = channel_labels(chani);  
% channel data for 1 trial
data    = EEG.data(chani,:,trial2plot);
% initialize Fourier output matrix
fourier = zeros(size(data)); 

%N= actual frequencies in Hz that will be returned by the
% Fourier transform. The number of unique frequencies we can measure is
% exactly 1/2 of the number of data points in the time series (plus DC). 

% Fourier transform is dot-product between sine wave and data at each frequency
disp(['Calculating Fourier coefficients for channel:' sensor2use ' for Frequency=1, to Frequency=' num2str(N)])
    
for fi=1:N
    disp(['Calcuating FT=sin*chani for freq=' num2str(fi) ' / ' num2str(N)]);
    sine_wave   = exp(-1i*2*pi*(fi-1).*time);
    fourier(fi) = sum(sine_wave.*data);
end

disp(['DONE with Fourier coefficients for channel:' sensor2use])
%normalize fourier
fourier=fourier/N;
% Save power spectrum for channel in mat file
% PowerSpectr_filename.mat
%subject-condition, channel, angle(fourier(1:N/2+1)),
%abs(fourier(1:N/2+1)), mean and std
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
powerspecfile = strcat('PowerSpectra_',mat_name); 
powerspecfile = strcat(powerspecfile,'.mat'); 

% abs fourier coeff
it_angleft= angle(fourier(1:N/2+1));
it_absft = abs(fourier(1:N/2+1)).^2;
absft(chani-1,1:end-1)= it_absft;
mean_absft(chani-1) = mean(it_absft);
std_absft(chani-1) = std(it_absft);
%angle fourier coeff
angleft(chani-1,1:end-1) = it_angleft;
mean_angleft(chani-1) = mean(it_angleft);
std_angleft(chani-1) = std(it_angleft);

%displaypowerspec = 1 to display power spectrum charts, one per channel, if
%0 only savemat file with data
displaypowerspec = 1;
if displaypowerspec == 1
    disp('Displaying the results')
    h = figure
    subplot(221)
    plot(data,'-o')
    set(gca,'xlim',[0 N+1])
    xlabel('Time')
    ylabel('Amplitude')
    title(['Time domain, Channel=' sensor2use]) 
    
    subplot(222)
    %plot3(frequencies_short,angle(fourier(1:N/2+1)),abs(fourier(1:N/2+1)).^2,'-o','linew',3)
    %shortend = 70*N/2; %only till fq = 70
    %freq_short = frequencies(1):frequencies(shortend);
    plot3(frequencies,it_angleft,it_absft,'-o','linew',3)
    grid on
    xlabel('Frequency (Hz)')
    ylabel('Phase')
    zlabel('Power')
    title('3-D of FT')
    view([20 20])
    
    subplot(223)
    bar(frequencies,it_absft)
    set(gca,'xlim',[-5 105])
    xlabel('Frequency (Hz)')
    ylabel('Power')
    title('Power spectrum discrete FT')
    
    subplot(224)
    bar(frequencies,it_angleft)
    set(gca,'xlim',[-5 105])
    xlabel('Frequency (Hz)')
    ylabel('Phase angle')
    set(gca,'ytick',-pi:pi/2:pi)
    title('Phase spectrum discrete FT')
    %savinf figure in figures directory
    disp('saving figure in figures directory')
    destfigfile = strcat('fig_power_channel_', sensor2use);
    destfigfile = strcat(destfigfile, '_');
    destfigfile = strcat(destfigfile, mat_name); 
    destfigfile = strcat(destfigfile,'.fig');  
    if ~exist('figures', 'dir')
         mkdir('figures');
    end
    olderfolder = cd('figures');
    savefig(h,destfigfile);
    cd(olderfolder); 
   
end

end
disp(['Saving mat file:' powerspecfile ' in folder:' mat_pathstr]); 
save(powerspecfile,'condsubjdate','absft','mean_absft','std_absft','angleft','mean_angleft','std_angleft')
disp('DONE!')   
end

function [] = remaincharts()
%% Figure 11.6

% Compute sine waves and sum to recover the original time series
reconstructed_data = zeros(size(data));
for fi=1:N
    % scale sine wave by fourier coefficient
    sine_wave = fourier(fi)*exp(1i*2*pi*(fi-1).*time);
    % sum sine waves together (take only real part)
    reconstructed_data = reconstructed_data + real(sine_wave);
end

figure
plot(data,'-o')
hold on
plot(reconstructed_data,'r-*')

legend({'original data';'inverse Fourier transform data'})

%% Figure 11.7

fft_data = fft(data)/N;

figure
subplot(131)
plot(frequencies,abs(fourier(1:N/2+1)).^2,'*-')
hold on
plot(frequencies,abs(fft_data(1:N/2+1)).^2,'ro-','markersize',8)
% make plot look nice
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power spectrum derived from discrete Fourier transform and from FFT')
axis square
legend({'time-domain Fourier';'FFT'})

subplot(132)
plot(frequencies,angle(fourier(1:N/2+1)),'*-')
hold on
plot(frequencies,angle(fft_data(1:N/2+1)),'ro-','markersize',8)
% make plot look nice
xlabel('Frequency (Hz)')
ylabel('Phase')
set(gca,'ytick',-pi:pi/2:pi)
title('Phase spectrum derived from discrete Fourier transform and from FFT')
axis square

subplot(133)
plot(reconstructed_data,'*-')
hold on
plot(ifft(fft(data)),'ro-','markersize',8)
% make plot look nice
xlabel('Time')
ylabel('Amplitude')
title('Manual inverse Fourier transform and ifft')
axis square

%% Figure 11.9

% list some frequencies
frex = [ 3 10 5 7 ];

% list some random amplitudes
amplit = [ 5 15 10 5 ];

% phases... 
phases = [  pi/7  pi/8  pi  pi/2 ];

% create a time series of sequenced sine waves
srate = 500;
time = -1:1/srate:1;
stationary  = zeros(1,length(time)*length(frex));
nonstationary = zeros(1,length(time)*length(frex));

for fi=1:length(frex)
    
    % compute sine wave
    temp_sine_wave = amplit(fi) * sin(2*pi*frex(fi).*time + phases(fi));
    
    % enter into stationary time series
    stationary = stationary + repmat(temp_sine_wave,1,length(frex));
    
    % optional change of amplitude over time
    temp_sine_wave = temp_sine_wave.*(time+1);
    
    % determine start and stop indices for insertion of sine wave
    start_idx = (fi-1)*length(time)+1;
    stop_idx  = (fi-1)*length(time)+length(time);
    
    % enter into non-stationary time series
    nonstationary(start_idx:stop_idx) = temp_sine_wave;
end

figure

% plot stationary signal
subplot(221)
plot(stationary,'r')
set(gca,'xlim',[1 length(stationary)])
title('stationary signal')

% plot non-stationary signal
subplot(222)
plot(nonstationary)
set(gca,'xlim',[1 length(nonstationary)])
title('non-stationary signal')

% perform FFT and plot
frequencies       = linspace(0,srate/2,length(nonstationary)/2+1);
fft_nonstationary = fft(nonstationary)/length(nonstationary);
fft_stationary    = fft(stationary)/length(stationary);
subplot(212)
plot(frequencies,abs(fft_stationary(1:length(frequencies)))*2,'r')
hold on
plot(frequencies,abs(fft_nonstationary(1:length(frequencies)))*2)
set(gca,'xlim',[0 max(frex)*2])
legend({'Power stationary';'Power non-stationary'})

%% Figure 11.10

% these figures produce the unassembled components of figure 10

load sampleEEGdata
eegdat4convol = squeeze(EEG.data(47,:,1));

% create Gaussian (you'll learn more about this formula in the next chapter)
time = -1:1/EEG.srate:1;
s = 5/(2*pi*30);
gaussian = exp((-time.^2)/(2*s^2))/30;

figure
subplot(211)
plot(eegdat4convol)


subplot(212)
plot(gaussian)


figure
subplot(211)
plot(conv(eegdat4convol,gaussian,'same'))

subplot(212)
plot(abs(fft(conv(eegdat4convol,gaussian,'same'))))

figure
subplot(211)
plot(abs(fft(eegdat4convol)))

subplot(212)
plot(abs(fft(gaussian)))

%% Figure 11.11

srate = 1000;
time  = -.5:1/srate:.5-1/srate;
f     = 20;
fg    = [15 5];
s     = sin(2*pi*f*time);

for i=1:2
    
    % compute Gaussian
    g = exp((-time.^2)/(2*(4/(2*pi*fg(i))^2)))/fg(i);
    
    
    figure
    
    subplot(411)
    plot(time,s)
    title('Sine wave (signal)')
    set(gca,'ylim',[-1.1 1.1])
    
    subplot(412)
    plot(time,g)
    title('Gaussian (kernel)')
    
    subplot(413)
    plot(time,conv(s,g,'same'))
    set(gca,'ylim',[-1.1 1.1])
    title('result of convolution')
    
    subplot(427)
    fft_s = abs(fft(s));
    fft_s = fft_s(1:floor(length(fft_s)/2)+1)./max(fft_s(1:floor(length(fft_s)/2)+1));
    bar(0:500,fft_s,'r')
    hold on
    
    fft_g = abs(fft(g));
    fft_g = fft_g(1:floor(length(fft_g)/2)+1)./max(fft_g(1:floor(length(fft_g)/2)+1));
    plot(0:500,fft_g)
    set(gca,'xlim',[0 40],'ylim',[0 1.05])
    title('individual power spectra')
    
    subplot(428)
    bar(0:500,fft_g.*fft_s)
    set(gca,'xlim',[0 40],'ylim',[0 .035])
    title('multiplied power spectra')
end

% inset scaling: axis([15 25 -.01 .11])

%% Figure 11.12

% create Gaussian
time = -1:1/EEG.srate:1;
s = 5/(2*pi*30);
gaussian = exp((-time.^2)/(2*s^2))/30;

figure

% plot EEG data
subplot(411)
plot(EEG.times,eegdat4convol)

% plot Gaussian
subplot(412)
plot(time,gaussian)

% plot result of convolution
subplot(413)
plot(EEG.times,eegdat4convol,'r')
hold on
plot(EEG.times,conv(eegdat4convol,gaussian,'same'))

subplot(427)
nfft = length(eegdat4convol);
fft_s = abs(fft(eegdat4convol,nfft));
fft_s = fft_s(1:floor(nfft/2)+1);
f = linspace(0,EEG.srate/2,floor(nfft/2)+1);
plot(f,fft_s./max(fft_s),'r')
hold on

fft_g = abs(fft(gaussian,nfft));
fft_g = fft_g(1:floor(nfft/2)+1);
plot(f,fft_g./max(fft_g))

set(gca,'xlim',[0 60])

subplot(428)
plot(f,fft_s.*fft_g)
set(gca,'xlim',[0 60])

%% end.
end