function [myfullname, EEG_study, channel_labels, initchan, centerfreq, chantouse, time] =  initialize_EEG_variables()
clear all
mydir = 'C:\Users\Jaime\Documents\BIAL PROJECT\patients\twh27\'
myfile = 'EEG_cut_BL_HYP_bs27_10222015_s2.mat'
myfullname = fullfile(mydir, myfile);
disp('Loading mat file...')
load(myfullname)
disp(['File ' myfile ' loaded!' ])

%rename the EEG object
EEG_study = EEG_cut_BL_HYP;
channel_labels = zeros(EEG_study.nbchan);
channel_labels = {EEG_study.chanlocs.labels};
disp([ 'Displaying the label of all the channels....' ])
initchan = 2 % channel 1 is the null Event channel
for i=initchan:EEG_study.nbchan
    disp(['Index ' i ' is the channel' channel_labels(i)])
end
chantouse = EEG_study.nbchan -1;
centerfreq = 10; % in Hz
%time = -1:1/EEG.srate:1;
time = EEG_study.times;
%initialize wavelet
end

