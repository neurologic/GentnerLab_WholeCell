function [upinds, lowinds] = WindowResponse(response_data, ci, windowsize, binsize, p)
%works best if windowsize is integer multiple of binsize
%this function takes a vector of data from which windows 
% of continuous to mostly continuous values above or below some confidence interval 
%the window (windowsize in samples) slides over the data sample by sample
%the window is binned (binsize in samples) and the mean value in each bin is calculated
%the window is significant if a critical number of bins within the window
%at that time point meet criteria
%p is the proportion of bins in window that must meet criteria
if size(response_data,2) < size(response_data,1)
    error('response_data must be column data')
end

win = [0,windowsize-1];
bins = [1:binsize:windowsize];

upinds = [];
lowinds = [];
for isamp = 1:size(response_data,2)- windowsize
    thiswin = win+isamp;
    windata = response_data(1,thiswin(1):thiswin(2));
    binvalue = [];
    for ibin = 1:(size(bins,2))
       binvalue(ibin) = mean(windata(1,bins(ibin):(bins(ibin)+binsize-1))); 
    end
    p_up = size(find(binvalue > ci(2)),2) / size(bins,2);
    p_low = size(find(binvalue < ci(1)),2) / size(bins,2);
    if p_up > p
        upinds = union(upinds, [thiswin(1):thiswin(2)]);
    end
    if p_low > p
        lowinds = union(lowinds, [thiswin(1):thiswin(2)]);
    end
end

%for now just want to return indices, not window edges, because want to
%make edges across db, not just for this signal/response
%
%
% d_up = diff(upinds)';
% break_up = find(d_up > 2);
% windows=min(upinds);
% for ibrk=1:size(break_up,2)
%     windows=[windows,upinds(break_up(ibrk)),upinds(break_up(ibrk)+1)];
% end
% windows=[windows,max(upinds)];
% windows=reshape(windows',2,size(windows,2)/2);
% win_up=windows;
% 
% d_low = diff(lowinds)';
% break_low = find(d_low > 2);
% windows=min(lowinds);
% for ibrk=1:size(break_low,2)
%     windows=[windows,lowinds(break_low(ibrk)),lowinds(break_low(ibrk)+1)];
% end
% windows=[windows,max(lowinds)];
% windows=reshape(windows',2,size(windows,2)/2);
% win_low=windows;
% 
% for iresp=1:size(windows,2)
% SigTimeBox(gca, windows(1,iresp),windows(2,iresp), get(gca,'YLim'),'r');
% end
% SigTimeBox(gca,1,size(response_data,2), [ci(1),  ci(2)],'k');
