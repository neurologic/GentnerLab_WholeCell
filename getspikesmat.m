function [spikesmat, gausstosmooth]=getspikesmat(highpassdata,threshold,dt)
spikesmat=[];

a=highpassdata>=threshold;
adif=diff(a')';
spikesmat=adif==1;

if isempty(spikesmat)
    return
end


for itrial=1:size(spikesmat,1)
    spiketimes{itrial}=find(spikesmat(itrial,:))*dt;
    meanisi(itrial)=mean(diff(spiketimes{itrial}))/dt;
    
end

allisi=round(nanmean(meanisi)/2);
if isnan(allisi)
%     error='no spikes in any trial'
%     return
gausstosmooth=[];
end
gausstosmooth=[];
if ~isnan(allisi);
    %     gausstosmooth = fspecial('gaussian', [allisi*4,1],allisi/2);
    gausstosmooth = fspecial('gaussian', [allisi*4,1],allisi/4);
    gausstosmooth = gausstosmooth/max(gausstosmooth);
end
