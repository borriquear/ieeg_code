%% Script for transforming the EDF from Natus acquisition system
%%into .mat file 
% JLPV Lab, Toronto ON 2015
% Gomez Ramirez, Jaime. PhD
%
%% Load the EDF of the session for Analysis
% pop_biosig()
% choose channel list, separated with commas e.g. [1,2,3,67,70]
%Invoke the functions
%n = number of channels
%times in function time_ini, time_end edftomat
%[list, t1,t2] = edftomat(n)
function [listofchannels, seconds_ini, seconds_end] = edftomat(n)
    %transform n in a list of numbers seprated by comas to be paste in
    %pop_biosig() to specify the list of channels to load in memory
    listofchannels = enumeratechannels(n);
    %print listlistofchanelslistofchannelsofchanels separated by commas
    allOneString = sprintf('%.0f,' , listofchannels);
    allOneString = allOneString(1:end-1);% strip final comma
    sprintf(allOneString)
    time_ini = '13:55:08'; %'12:53:05 ';  induction timelimits = [50108 50264] || normal baseline[52559 52799]
    time_end = '13:57:44'; %12:56:41
    [seconds_ini] = calcseconds(time_ini);
    [seconds_end] = calcseconds(time_end);
    epoch_seconds = seconds_end - seconds_ini;
    sprintf('Seconds_init = %.0f, Seconds_end = %.0f, Seconds_divergence = %.0f',seconds_ini,seconds_end,epoch_seconds)
end
function lofch = enumeratechannels(n)
    lofch = [];
    for i=1:n
        lofch = [lofch i];
    end  
end

function sec_time = calcseconds(day_time)
    [Y, M, D, H, MN, S] = datevec(day_time);
    sec_time = H*3600+MN*60+S;
end
