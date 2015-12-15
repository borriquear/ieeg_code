%% Script for preprocesssing iEEG data from Natus acquisition system
% JLPV Lab, Toronto ON 2015
% Gomez Ramirez, Jaime. PhD
%
%% Pre processing of mat file 
clear all
patinetsfolder = 'C:\Users\Jaime\Documents\BIAL PROJECT\patients\'
singlepatientfolder = 'Cordeiro, Jessica'
sessionspfolder = 'Session 1 day 1 20 October'
sessionmatfile = 'test-10302015-all-10minsr250.mat'
lfile = fullfile(patinetsfolder,singlepatientfolder,sessionspfolder,sessionmatfile)
if exist(lfile, 'file') == 2
    load(lfile)
else
    % File does not exist.
    warningMessage = sprintf('Warning: file does not exist:\n%s', lfile);
    uiwait(msgbox(warningMessage));
    return;
end
%
% Optionally, open dialog to open file this would also work:
% [file2load,path4file]=uigetfile('*.mat','Please select EEG data file');
% load([ path4file file2load ])

%The EEG structure is loaded. Check number of channels and other fields
%EEG.nbchan
%plot a couple of channels
%plot(EEG.times, EEG.data(1,:,:))
%plot(EEG.times,EEG.data(2,:,:))

%% Filtering the signal
% resample the signal is necessary
%pop_resample()
locutoff = 0.5 %high band filter
hicutoff = 70  %low band filter
notchf = 60 %notch filter in the Americas 
%http://ftp.aas.duke.edu/pub/pc/classrooms/eeglab12_0_2_3b/plugins/firfilt1.5.4/pop_eegfiltnew.m
sprintf('Filtering signal, High pass =%0.1f Hz. and Low Band %0.1f = Hz.',locutoff,hicutoff )
[EEG, com, b] = pop_eegfiltnew(EEG, locutoff, hicutoff)
% notch filter at 60 Hz.the filter passes all frequencies, except for the range of 59–61 Hz
sprintf('Filtering with notch filter at %0.1f Hz.',notchf )
EEG = pop_iirfilt( EEG, 59, 61, [], [1]);
% Save EEG data into .mat file
filteredsessionmatfile = 'fil_test-10302015-all-10minsr250.mat'
save(filteredsessionmatfile,'EEG')

% Plot channel spectra maps
%% Removing bad channels
pop_eegplot(EEG)

%% Removing bad epochs


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

%%

