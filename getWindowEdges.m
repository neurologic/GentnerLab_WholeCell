function windows = getWindowEdges (includeinds, timebetween, timewithin)
%this function takes a vector of discontinuous indices and returns the
%windows where the indices are continuous
%the windows are constrained by lumping any windows that are less than
%"timebetween" apart and keeping only windows that are "timewithin" long
%timebetween and timewithin are in samples (indices)
if size(includeinds,2) < size(includeinds,1)
   includeinds = includeinds'; 
end
tmp=sort(includeinds);
d=diff(tmp);


respbreakind=find(d>timebetween);
windows=min(includeinds);
for ibrk=1:size(respbreakind,2)
    windows=[windows,includeinds(respbreakind(ibrk)),includeinds(respbreakind(ibrk)+1)];
end
windows=[windows,max(includeinds)];
windows=reshape(windows',2,size(windows,2)/2);

resptimes=windows(2,:)-windows(1,:);
%each response must last longer than 50ms to count as a response
%window
useresp=find(resptimes>timewithin);
windows=windows(:,useresp);