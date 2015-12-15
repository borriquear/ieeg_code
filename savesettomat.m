% convert .set into .mat

 mydir = 'C:\Users\shelagh\Desktop\patients\Cordeiro, Jessica\Session 1 day 1 20 October\dp_cj28_20151020_s1'
 myfile = 'EEG_raw_cj28_20151020_s1.mat'
 myfullname = fullfile(mydir, myfile)
 save(myfullname, 'EEG')
 timelimits2 = [52799 5526122]
 EEGcutbsnormal =  eeg_eegrej(EEGcutbsnormal1, timelimits2);
 save('cut1_bs_induction_filteredbp-notch_cj28_20151020_s1.mat','EEGcutbsinduction')