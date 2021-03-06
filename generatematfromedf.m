%% Preprocesssing iEEG data from Natus acquisition system
% Author: Jaime Gomez Ramirez, Ph.D., mathematical neuroscience
% The Hospital for Sick Chldren, Toronto, ON
% email address: jd.gomezramirez@gmail.com  
% Website: aslab.org/~gomez
% November 2015; Last revision: 16-February-2015
%
function [listofchannels] = generatematfromedf(n)
%% edftomat given the edf file obtains the mat file filtered .
% NOTE this function is not done, and this process is doneusing the GUI of
% eeglab
% cuteonepoch function assumes that the mat is created and cut one epoch
% for the mat file of the entore session
% So far this functio only does very basic stuff like generate the list of
% channels

%transform n in a list of numbers seprated by comas to be paste in
%pop_biosig() to specify the list of channels to load in memory
listofchannels = enumeratechannels(n);

preprocessingofedf();
%print listlistofchanelslistofchannelsofchanels separated by commas
% allOneString = sprintf('%.0f,' , listofchannels);
% allOneString = allOneString(1:end-1);% strip final comma
% sprintf(allOneString)
end
function preprocessingofedf()
%% This function needs to be check.
% Tries to do the fultering automatically

%% Load the edf file 
patinetsfolder = 'C:\Users\shelagh\Desktop\patients'
singlepatientfolder = 'Cordeiro, Jessica'
sessionspfolder = 'Session 1 day 1 20 October'
dpfolder = 'dp_cj28_20151020_s1'
sessionmatfile = 'EEG_raw_cj28_20151020_s1.mat'
lfile = fullfile(patinetsfolder,singlepatientfolder,sessionspfolder, dpfolder, sessionmatfile)
if exist(lfile, 'file') == 2
    load(lfile)
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', lfile);
    uiwait(msgbox(warningMessage));
    return;
end
%% Filtering and resampling
% Optionally, open dialog to open file this would also work:
% [file2load,path4file]=uigetfile('*.mat','Please select EEG data file');
% load([ path4file file2load ])
% Resampling and  Filtering the signal
% Resample the signal is necessary
% Optionally we can resample the EEG 
newfreq = 256
[EEG] = pop_resample(EEG, newfreq);
filtered_rsample_sessionmatfile = 'EEG_raw_rsample256_cj28_20151020_s1.mat'
lfiltfile = fullfile(patinetsfolder,singlepatientfolder,sessionspfolder, dpfolder,filtered_rsample_sessionmatfile);
save(lfiltfile,'EEG')
%Filter data band pass and stop filters
locutoff = 0.5 %high band filter
hicutoff = 70  %low band filter
notchf = 60 %notch filter in the Americas 
%http://ftp.aas.duke.edu/pub/pc/classrooms/eeglab12_0_2_3b/plugins/firfilt1.5.4/pop_eegfiltnew.m
sprintf('Filtering signal, High pass =%0.1f Hz. and Low Band %0.1f = Hz.',locutoff,hicutoff )
[EEG, com, b] = pop_eegfiltnew(EEG, locutoff, hicutoff)
% notch filter at 60 Hz.the filter passes all frequencies, except for the range of 59�61 Hz
sprintf('Filtering with notch filter at %0.1f Hz.',notchf )
EEG = pop_iirfilt( EEG, 59, 61, [], [1]);
% Save EEG data into .mat file
filteredsessionmatfile = 'EEG_raw_filtered_rsample256_cj28_20151020_s1.mat'
lfiltfile = fullfile(patinetsfolder,singlepatientfolder,sessionspfolder, dpfolder,filteredsessionmatfile);
save(lfiltfile,'EEG','-v7.3')

% Plot channel spectra maps
%% Removing bad channels
%pop_eegplot(EEG)

%% Removing bad epochs

%% Hilbert Transform
%load sampleEEGdata
% ERP, and its hilbert transform
erp  = squeeze(mean(EEG.data(2,:,:),3));
erpH = hilbert(erp);

figure(2), clf

% plot ERP and real part of Hilbert transformed ERP
subplot(311)
plot(EEG.times,erp), hold on
plot(EEG.times,real(erpH),'r')
legend({'ERP';'real(hilbert(erp))'})
legend({'real(hilbert(erp))'})

% plot ERP and magnitude
subplot(312)
plot(EEG.times,erp), hold on
plot(EEG.times,abs(erpH),'r')
legend({'real(hilbert(erp))';'abs(hilbert(erp))'})

% plot ERP and phase angle time series
subplot(313)
plot(EEG.times,erp), hold on
plot(EEG.times,angle(erpH),'r')
legend({'real(hilbert(erp))';'angle(hilbert(erp))'})
xlabel('Time (ms)'), ylabel('Voltage or radians')

% plot as 3d line
figure(3), clf
plot3(EEG.times,real(erpH),imag(erpH))
xlabel('Time (ms)'), ylabel('Real part'), zlabel('Imaginary part')
axis tight
rotate3d

%% finding time indices based on ms

% We want to create a topographical plot at time=300 ms
% time index = sratexseconds, e.g. to plot the topography at 27 seconds 
%
time2plot = 300; % in ms!
[minval, minidx] = min(abs(EEG.times - time2plot))
% extract the trial-averaged data from requested time point
data2plot = squeeze(mean( EEG.data(:,minidx,:) ,3));

% plot
figure(1), clf
topoplot(data2plot,EEG.chanlocs);
title([ 'Topoplot from time=' num2str(EEG.times(time2plot)) ' ms.' ])

%% same concept for frequencies
% vector of 42 linearly spaces freqs (in Hz) betw 2 and 100
frex = linspace(2,100,42);
%if we want to analyzedata at 23 Hz.
freqIwant = 23; % in hz

% use min(abs trick to find closest frequency to 23 Hz
[~,frexidx] = min(abs(frex-freqIwant));

% the function dsearchn also works
frexidx = dsearchn(frex',freqIwant);

%% indexing channels based on names
% to see all the channeks labels : >> {EEG.chanlocs.labels}
%EEG.chanlocs(2): 
%         labels: 'Grid1'
%            ref: ''
%          theta: []
%         radius: []
%              X: []
%              Y: []
%              Z: []
%      sph_theta: []
%        sph_phi: []
%     sph_radius: []
%           type: ''
%         urchan: []
% the electrode label that we want to analyze
electrodeName = 'p1'; % case doesn't matter

% find the channel number that corresponds to this label
electrodeidx = strcmpi(electrodeName,{EEG.chanlocs.labels});

% confirm that this electrode is correct
EEG.chanlocs(electrodeidx)

% plot the ERP from this electrode
figure(1), clf
plot(EEG.times,mean( EEG.data(electrodeidx,:,:),3 ))

%% indexing multiple electrodes

electrodeNames = {'Grid1';'Grid2';'Grid3'};
% initialize
electrodeidx = zeros(size(electrodeNames));

% loop through electrodes and find the index of each one
% (The code below has three errors. Try to find/fix them!)
for chani=1:length(electrodeNames)
    try 
        electrodeidx(chani) = find(strcmpi(electrodeNames{chani},{EEG.chanlocs.labels}));
    catch me;
        warning('Error with the labeling of the electodes');
    end
end
% plot all the ERPs
%for onlyone trial plot(EEG.times, EEG.data(electrodeidx,:,:))
plot(EEG.times,mean( EEG.data(electrodeidx,:,:),3 ))
legend({EEG.chanlocs(electrodeidx).labels})
end


function lofch = enumeratechannels(n)
%% enumeratechannels create a list of numbers 1..n 
% in: n number of channels
% output lofch 
% paste lofch within [] when loading the edf file for those channel numbers

    lofch = [];
    initchan =1;
    for i=initchan:n
        lofch(i) =  i;
    end  
end
