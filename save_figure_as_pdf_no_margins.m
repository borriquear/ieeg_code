function [  ] = save_figure_as_pdf_no_margins( filelist)
%%save_figutr_as_pdf_no_margins 
% if filelist = [] save for specific figure index that is open
% else filelist contains the list of the fig fies (entire path) to be saved as
% pdf

if filelist == []
    
    figureindex =1 % number pf the figure we want to save
    h = figure(figureindex);
    ti = get(gca,'TightInset')
    set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);
    set(gca,'units','centimeters')
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
    saveas(h,'output.pdf');
    %% Save to pdf from .fig file in computer
    path_fig_files = 'C:\Users\Jaime\Documents\BIAL PROJECT\patients\figures';
    filelist = getAllFiles(path_fig_files);
else

end

end


