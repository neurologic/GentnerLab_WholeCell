saveme=0;
allVm=[999];
numrepstims=1;
havestimcond=expt.stimcond;

plotrasterinfo.isplotraster=1;
plotrasterinfo.spikesVm=999;
plotrasterinfo.spikesthresh=10;
plotrasterinfo.cutoff=100;
plotrasterinfo.scaleticks=5;
plotrasterinfo.negative=1;

plotrestricted.isplotrestricted=1;
plotrestricted.AfterMotDur=1.5;

lowpass=1;
meanonly=1;
ncol=2;

for iVm=1:length(allVm)
    subplotMotifByVm(expt,allVm(iVm),ncol,lowpass,params,havestimcond,numrepstims,r.Dir.Wavs,plotrasterinfo,plotrestricted,'mV',r,saveme,meanonly)
end

%% filtering out trials when had the juxtacellular recording
cutoff=100;
highpassdata=HighpassGeneral(expt.wc.data,expt.analysis.params.baselinewin,[],cutoff,1/expt.wc.dt);

threshold=75;
bintime=0.05;
spikesmat=getspikesmat(bintime,highpassdata,expt.wc.dt,threshold);

ntrials=size(highpassdata,1);

triallength=size(highpassdata,2)*expt.wc.dt;
index=1;
for itrial=1:ntrials
    spiketimes=find(spikesmat(itrial,:));
    if (size(spiketimes,2)/triallength)>10
        trialinds(index)=itrial;
        index=index+1;
%     else
%         trialinds(itrial)=0;
    end
end

Lexpt=filtesweeps(expt,0,'trial',trialinds);