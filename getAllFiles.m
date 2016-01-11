function fileList = getAllFiles(dirName)
%getAllFiles as pdf from directory
% get the list of files in fileList
dirData = dir(dirName);      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
        fileList,'UniformOutput',false);
end
subDirs = {dirData(dirIndex).name};  %# Get a list of the subdirectories
validIndex = ~ismember(subDirs,{'.','..'});  %# Find index of subdirectories
%#   that are not '.' or '..'
for iDir = find(validIndex)                  %# Loop over valid subdirectories
    nextDir = fullfile(dirName,subDirs{iDir});    %# Get the subdirectory path
    fileList = [fileList; getAllFiles(nextDir)];  %# Recursively call getAllFiles
end
% convert get the list of files in fileList to char, and save each fig as
% pdf
%
for i=1:size(fileList)
    hfigfile = char(fileList(i));  
    [pathstr,name,ext] = fileparts(hfigfile); 
    dest_pdf_file = strcat(name, '.pdf');
    full_dest_file = fullfile(pathstr,dest_pdf_file);
    disp([' opening...' num2str(i)])
    hfig = openfig(hfigfile, 'invisible');
    ti = get(gca,'TightInset');
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    set(gca,'units','centimeters');
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    saveas(hfig,full_dest_file);
    %close(hfig);
end
disp('closing all')
end