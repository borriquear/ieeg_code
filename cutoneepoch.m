%% Preprocesssing iEEG data from Natus acquisition system
% Author: Jaime Gomez Ramirez, Ph.D., mathematical neuroscience
% The Hospital for Sick Chldren, Toronto, ON
% email address: jd.gomezramirez@gmail.com  
% Website: aslab.org/~gomez
% November 2015; Last revision: 16-February-2015

function [] = cutoneepoch()
%% preprocessing generate one epoch segment from a mat file of the entire session 
% DEPENDENCIES: It needs the .mat file of the entire session filtered in
% localdir/localmat
% the .mat file has to be already been generated from the corresponding EDF file and
% and band and notch filtered 

%% Load the .mat file with the entire session. This section can be skipped if we load the file form the GUI
localdir = 'D:\BIAL PROJECT\patients\On, Freeman'
localmat = 'dummy1chan-filternadnotch.mat'
%patient , date and session must correspond with the patient of the entire sessiondefined in localmat
patientlist = {'CJ28' 'CJ28' 'TWH29' 'TWH30' 'TWH31' 'TWH32' 'TWH33' 'TWH34' 'TWH35'};
patientindex = 1;
se_patient = patientlist(patientindex); 
se_date = '20151020'; % MMDDYYYY
nu_session = 's1'; %sN
fullmatfile = fullfile(localdir, localmat);
disp(['Loading mat file alrady filtered ...' fullmatfile])
% load the .mat file that contains the entire sesssion already filtered
load(fullmatfile);
EEGentire = EEG;

%% Specify the times to cut inthe session file initially load
% check excel file to see the initial and final time of the epoch

hmsdate_startsession = '12:26:07';
hmsdate_initepoch = '12:53:00' ;
hmsdate_endepoch = '12:54:42';
% get the seconds from the actualtime
% t0...t10.....t20, t21 is the duration in seconds 
EEG = EEGentire;
disp(sprintf('From time to seconds, init session= %s , initial epoch= %s respectively', hmsdate_startsession, hmsdate_initepoch ))
[t10] = fromtime_to_seconds(hmsdate_startsession, hmsdate_initepoch);
disp(['seconds, t10 = ', num2str(t10)])
[t21] = fromtime_to_seconds(hmsdate_initepoch, hmsdate_endepoch);
disp(['seconds, t21 = ', num2str(t21)])
t20 = t21 + t10;
timelimits = [t10 t20];
disp(['seconds t20=', num2str(t20)])
%scale timelimits with the samploing rate 
datapointtimelimits = timelimits*EEG.srate;

disp(['Creating epoch between [t10, t20]s=' num2str(timelimits)])
EEGepoch =  eeg_eegrej(EEGentire, datapointtimelimits);
disp(['EEGepoch created '])
% save .mat file
disp(['Saving EEGepoch as a mat file '])
savematfileforepoch(localdir, EEGepoch, se_patient, se_date, nu_session)
disp(['Created mat file for :' se_patient, se_date, nu_session])
end

function [] = savematfileforepoch(localdir, EEGepoch, se_patient, se_date, nu_session)
%% savematfile save mat file for the patient, condition etc specified in the function
% INPUT: directory where mat file is saved
% EEG object
% se_patient, se_date, nu_session to write the destination mat file name
%disp(['Patient, date and session are:' num2str(se_patient) , num2str(se_date), num2str(nu_session) ])
EEG = EEGepoch
matfh = 'EEG_cut_';
conditionlist = {'BL_EC_POST' 'BL_EO_POST' 'BL_EC_PRE' 'BL_EO_PRE' 'BL_HYP'};
conditionindex = 5;
%patientlist = {'TWH28' 'TWH29' 'TWH30' 'TWH31' 'TWH32' 'TWH33' 'TWH34' 'TWH35'  };
%patientindex = 1;
%se_date = '10202015'; % MMDDYYYY
%nu_session = 's1'; %sN
disp(['Condition is:' conditionlist(conditionindex)])
% matfh = strcat(matfh, conditionlist(conditionindex));
% matfh = strcat(matfh, '_');
% matfh = strcat(matfh, se_patient);
% matfh = strcat(matfh, '_');
% matfh = strcat(matfh, se_date);
% matfh = strcat(matfh, '_');
% matfh = strcat(matfh, nu_session);
matfh = sprintf('EEG_cut_%s_%s_%s_%s', conditionlist{conditionindex},se_patient{1}, se_date, nu_session)
disp(['mat file is:' matfh])
fullmatfile = fullfile(localdir, matfh);
oldFolder = cd(localdir)
disp(['Saving EEG epoch in file:', matfh])
save(matfh,'EEG') 
end

function [total_seconds_difference] = fromtime_to_seconds(hmsdate_first,hmsdate_last)
%FUNCTION_NAME - function that takes two strings with format hour:minute:second and returns the seconds, the first is the initial time and the second is the last time 
%and calculates the number of seconds that passed between the two instants 
%
% Inputs:
%    hmsdate_first - literal with the form 'h:m:s' for initial time
%    hmsdate_last -  literal with the form 'h:m:s' for final time
% Outputs:
%    total_seconds_difference  = sesoncdsinit- seconds final
% Example: 
%   [output1] = fromtime_to_seconds('12:34:59', '12:39:52')
%
[Y, M, D, H, MN, S] = datevec(hmsdate_first);
total_seconds_init = H*3600+MN*60+S;
[Y, M, D, H, MN, S] = datevec(hmsdate_last);
total_seconds_final = H*3600+MN*60+S;
total_seconds_difference = total_seconds_final - total_seconds_init;
end