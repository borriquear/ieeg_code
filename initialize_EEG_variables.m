function [myfullname, EEG_study, channel_labels] =  initialize_EEG_variables()
% initialize_EEG_variables  load mat file and initialize EEG lab variables .
%   [myfullname, EEG_study, channel_labels] =
%   initialize_EEG_variables() load matfile, save the EEG object as
%   EEG_study. Note this is hard-code: EEG_study = EEG_cut_BL_HYP
%   channel_labels labels of the channels
%   myfullname contains information of the condition,session, patient and date 
%

clear all
%Load matfile 
mydir = 'D:\BIAL PROJECT\patients\Way, SA\'
%mydir = 'D:\BIAL PROJECT\patients\twh31\'
myfile = 'EEG_cut_BL_HYP_was30_11172015_s1.mat'
%myfile ='EEG_cut_BL_HYP_sm31_12012015_s1.mat'
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

%time = -1:1/EEG_study.srate:1;
end

