% get [IDX,X] = KmeansPredict(experiment's responses and stimuli)

binsizeVec = [500:500:4000]; %samples

exptnames = {
    'KP_B130_131120_p1c2c.mat'};  %at first, just working on the first experiment from this cell...
%so only dealing with one stimulus at
%first... then can go through and add
%data from other stimuli recorded
%     'KP_B130_131120_p1c2b.mat'
%     'KP_B130_131120_p1c2c.mat'};
numbinind = 1;

numUniqueTemplate = [];
numUnique = [];
numBins = [];

hfig = figure;
hold on
% templatefig = figure;
% hold on
for ib = 1:size(binsizeVec,2)
    binsize = binsizeVec(ib);
    IDXhat = [];
    Chat = [];
    X = [];
    Xavg = [];
    xind = 1;
    for iexpt = 1:size(exptnames,1)
        thisexpt = exptnames{iexpt};
        load([r.Dir.Expt thisexpt]); %load the experiment .mat file
        vmexpt = filtesweeps(expt,0,'Vm',0); %filter out all trials that were not recorded under current clamp at 0 holding
        stimcond = expt.stimcond; %the structure with information about the stimuli played
        %require a certain number reps
    table=getClampTab(expt,{'clamp',0});
    keepsigs=reprequire(table,trials);
    thiscond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
    stimcond = thiscond;
    %     stimcond = stimcond(keepsigs);
    loudstim = [];
    not_db_ws = [];
    notwarpstim = [];
    for istim = 1:size(stimcond,2) %grab only loudest db of db stims... or non-db stims from this set of stimuli
        %         if isempty(regexp(stimcond(istim).wavnames,'ws'))
        %             notwarpstim = [notwarpstim,istim];
        %         end
        if isempty(regexp(stimcond(istim).wavnames,'d')) && isempty(regexp(stimcond(istim).wavnames,'ws'))
            not_db_ws = [notdbstim, istim];
        end
        
        if ~isempty(regexp(stimcond(istim).wavnames,'d'))
            if ~isempty(regexp(stimcond(istim).wavnames,'d80')) || ~isempty(regexp(stimcond(istim).wavnames,'db80'))
                loudstim = [loudstim,istim];
            end
        end
    end
    
    
    stiminds = union(not_db_ws, loudstim); %left with only "unique" stimuli....
    %     stiminds = union(stiminds, notwarpstim);
    if isempty(stiminds)
        continue
    end
        
        %bin the response for that stimulus into substimuli of "binsize" long
        for istim = 1:max(size(stiminds))
            thisind = stiminds(istim);
            stimwav = stimcond(thisind).wavs; %this is the stimulus wavform in case you want to see it
            %(yes... it is a repeat of a motif, but that could end up
            %being something interesting to note... because the
            %response during the second repeat is usually different
            %from the first repeat even though the stimulus is exactly
            %the same... context/history effects
            stimexpt = filtesweeps(vmexpt,0,'wavnames',stimcond(thisind).wavnames); %filtering out all trials that is not that stimulus
            sigon = expt.analysis.params.baselinewin(2);
            sigoff = sigon+round(size(stimcond(thisind).wavs,1)/44100/expt.wc.dt);
            sigdata = medfilt1(stimexpt.wc.data,200,[],2); %all trials in response to that stimulus
            sigdata = sigdata(:,sigon:sigoff);
            numbins = floor((sigoff-sigon)/binsize); %how many subresponses of "binsize" width can there be
            bins = [1:binsize:numbins*binsize]; %make the edges of the subresponse bins
            for ibin = 1:numbins-1 %parse the response into subresponses
                X = [X;sigdata(:,bins(ibin):(bins(ibin+1)-1))]; %stack all of the trials of all of the subresponses on top of each other
                IDXhat = [IDXhat;repmat(xind,size(sigdata,1),1)]; %keep track of what bin each trial "belongs to"
                [idx,c] = kmeans(sigdata(:,bins(ibin):(bins(ibin+1)-1)),1); %get an estimate of the "real" cluster center
                Chat (ibin,:) = c; %keep track of the centroid for each set of trials for each bin
                Xavg = [Xavg;mean(sigdata(:,bins(ibin):(bins(ibin+1)-1)))]; %this is the mean Vm response for each bin
                xind = xind +1;
            end
        end
    end
    
    %I added this "template" part in as an idea for maybe how to constrain the
    %clustering a little so that it was more specific, but I am not sure that
    %it helps and ultimately it might be too non-agnostic to the trials data (X matrix clustering) to use it
%     
%     for ik = 1:numbins-1
%         replicates = 10; %how many times to run kmeans
%         for irep = 1:replicates
%             [IDXtemplate, Ctemplate,sumd] = kmeans(Xavg,numbins-1,'replicate',5); %calculate the centroid of the average Vm response for each subresponse
%             templateD(irep,ik) = sum(sumd);
%         end
%     end
%     %     figure(templatefig)
%     %     subplot(5,8,ib)
%     %     hold on
%     numKhyp = [];
%     numKhyp(1) = 1;
%     for ik = 1:numbins-1
%         %         scatter(repmat(ik,replicates,1),templateD(:,ik))
%         if ik>1
%             numKhyp(ik) = kstest2(templateD(:,ik-1),templateD(:,ik),'Tail','smaller','Alpha',0.01);
%         end
%     end
%     isFail = find(numKhyp==0);
%     if isempty(isFail)
%         numUniqueTemplate(ib) =  numbins-1;
%     end
%     if ~isempty(isFail)
%         numUniqueTemplate(ib) = min(isFail);
%     end
%     
    
    for ik = 1:numbins-1
        replicates = 10; %how many times to run kmeans
        for irep = 1:replicates
            [IDX, C, sumd,D] = kmeans(X,ik,'replicates',replicates);
            %use replicates = 2 so that does not break so much on empty
            %clusters... will automatically at least wrap to a second
            %replication
            %results of kmeans... IDX is the indices of what cluster it belongs to... and C is the centroid locations
            %Ctemplate was used to set the starting points for the centroids
            %for the analsysis...
            totalD(irep,ik) = sum(sumd);
        end
    end
    figure(hfig)
    subplot(5,8,ib)
    hold on
    numKhyp = [];
    numKhyp(1) = 1;
    for ik = 1:numbins-1
        scatter(repmat(ik,replicates,1),totalD(:,ik))
        if ik>1
            numKhyp(ik) = kstest2(totalD(:,ik-1),totalD(:,ik),'Tail','smaller','Alpha',0.01);
        end
    end
    isFail = find(numKhyp==0);
    if isempty(isFail)
        numUnique(ib) =  numbins-1;
    end
    if ~isempty(isFail)
        numUnique(ib) = min(isFail);
    end
    
    numBins(ib) = numbins-1;
    
    %     numbinind = numbinind + 1;
end

title(' ')
xlabel('k clusters','FontSize',20)
ylabel('total sum distances','FontSize',20)

figure;
scatter(binsizeVec*expt.wc.dt,(numUnique./numBins),50,'k','fill')
ylabel('number of k at which all clusters different / total bins','FontSize',14)
xlabel('seconds per bin','FontSize',14)

%% organize by cluster and test all combinations
getnoisefloor = 0;
binsizeVec = [1,100,400,600,800,1000:1000:10000];
% exptnames = {
%     'KP_B130_131120_p1c2c.mat'};
% numestimates = 10;
numsubsamp = 5;
trials = 5;
numrep = 100;
numbinind = 1;
MaxTotKeyacc = [];
MaxAvgKeyacc = [];
MaxStimacc = [];
exptunqinds = [];
exptstimnames = [];
domean = 1;
getbaseline = 0;
%%


for iexpt=1:size(repexpts,1)
    
    for iunq = 1:size(unq_expts,2) 
        if ~isempty(regexp(repexpts{iexpt},unq_expts{iunq}))
            exptunqinds(iexpt) = iunq;
        end
    end 
    
    thisexpt=repexpts{iexpt};
    
    load([r.Dir.Expt thisexpt]); %load the experiment .mat file
    vmexpt = filtesweeps(expt,0,'Vm',0); %filter out all trials that were not recorded under current clamp at 0 holding
    stimcond = expt.stimcond; %the structure with information about the stimuli played
    
    %require a certain number reps
    table=getClampTab(expt,{'clamp',0});
    keepsigs=reprequire(table,trials);
    thiscond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
    stimcond = thiscond;
    %     stimcond = stimcond(keepsigs);
    loudstim = [];
    not_db_ws = [];
    notwarpstim = [];
    for istim = 1:size(stimcond,2) %grab only loudest db of db stims... or non-db stims from this set of stimuli
        %         if isempty(regexp(stimcond(istim).wavnames,'ws'))
        %             notwarpstim = [notwarpstim,istim];
        %         end
        if isempty(regexp(stimcond(istim).wavnames,'d')) && isempty(regexp(stimcond(istim).wavnames,'ws'))
            not_db_ws = [notdbstim, istim];
        end
        
        if ~isempty(regexp(stimcond(istim).wavnames,'d'))
            if ~isempty(regexp(stimcond(istim).wavnames,'d80')) || ~isempty(regexp(stimcond(istim).wavnames,'db80'))
                loudstim = [loudstim,istim];
            end
        end
    end
    
    
    stiminds = union(not_db_ws, loudstim); %left with only "unique" stimuli....
    %     stiminds = union(stiminds, notwarpstim);
    if isempty(stiminds)
        continue
    end
    allnames = [];
    for istim = 1:max(size(stiminds))
        thisind = stiminds(istim);
        allnames{istim} = stimcond(thisind).wavnames;
    end
    exptstimnames{iexpt} = allnames;

    %     for iestimates = 1:numestimates
    for ib = 1:size(binsizeVec,2)
        %         numsubsamp = 5;
        binsize = binsizeVec(ib);
        IDXhat = [];
        Chat = [];
        X = [];
        XmeanV = [];
        Xavg = [];
        xind = 1;
        %bin the response for that stimulus into substimuli of "binsize" long
        for istim = 1:max(size(stiminds))
            thisind = stiminds(istim);
            stimwav = stimcond(thisind).wavs; %this is the stimulus wavform in case you want to see it
            %(yes... it is a repeat of a motif, but that could end up
            %being something interesting to note... because the
            %response during the second repeat is usually different
            %from the first repeat even though the stimulus is exactly
            %the same... context/history effects
            stimexpt = filtesweeps(vmexpt,0,'wavnames',stimcond(thisind).wavnames); %filtering out all trials that is not that stimulus
            sigon = expt.analysis.params.baselinewin(2);
            sigoff = sigon+round(size(stimcond(thisind).wavs,1)/44100/expt.wc.dt);
            sigdata = medfilt1(stimexpt.wc.data,200,[],2);
            if getbaseline == 0
                sigdata = sigdata(:,sigon:sigoff); %all trials in response to that stimulus
            end
            if getbaseline == 1
                sigdata = sigdata(:,expt.analysis.params.baselinewin(1):sigon);
            end
            
            
            if getnoisefloor == 1
                for itrial = 1:size(sigdata,1)
                    sigdata(itrial,:) = sigdata(itrial,randperm(size(sigdata,2)));
                end
            end
            
            ntrials = size(sigdata,1);
            if getbaseline == 0
            numbins = floor((sigoff-sigon)/binsize); %how many subresponses of "binsize" width can there be
            end
            if getbaseline == 1
            numbins = floor((sigon-expt.analysis.params.baselinewin(1))/binsize); %how many subresponses of "binsize" width can there be
            end
            bins = [1:binsize:numbins*binsize]; %make the edges of the subresponse bins
            for ibin = 1:numbins-1 %parse the response into subresponses
                if domean == 0;
                    X = [X;sigdata(:,bins(ibin):(bins(ibin+1)-1))]; %stack all of the trials of all of the subresponses on top of each other
                end
                if domean == 1;
                    X = [X; mean(sigdata(:,bins(ibin):(bins(ibin+1)-1)),2)];
                end
                IDXhat = [IDXhat;repmat(xind,size(sigdata,1),1)]; %keep track of what bin each trial "belongs to"
                [idx,c] = kmeans(sigdata(:,bins(ibin):(bins(ibin+1)-1)),1); %get an estimate of the "real" cluster center
                Chat (ibin,:) = c; %keep track of the centroid for each set of trials for each bin
                Xavg = [Xavg;mean(sigdata(:,bins(ibin):(bins(ibin+1)-1)))]; %this is the mean Vm response for each bin
                xind = xind +1;
            end
        end
        
        %         ik = numbins-1;
        
        replicates = 10; %how many times to run kmeans
        
        %subsample X
        
        if size(unique(IDXhat),1) < numsubsamp
            continue
        end
        
        hits = [];
        stimacc= [];
        keyacc = [];
        totalkeyacc = [];
        avgkeyacc = [];
        for iestimate = 1:numrep
            allstims = unique(IDXhat);
            subsamps = randsample(allstims,numsubsamp);
            subsampinds = [];
            for isub = 1:size(subsamps,1)
                subsampinds = [subsampinds;find(IDXhat == subsamps(isub))];
            end
            Xsub = X(subsampinds,:);
            IDXhatsub = IDXhat(subsampinds,:);
            [IDX, C, sumd,D] = kmeans(Xsub,numsubsamp,'replicate',replicates); %,'start',...
            %repmat(Ctemplate,[1,1,replicates]));
            
            allclusters = unique(IDX);
            allstims = subsamps;
            %             allstims = [1:numsubsamp];
            
            allKeys = perms(allclusters);
            
            for ikey = 1:size(allKeys,1)
                thiskey = allKeys(ikey,:);
                for istim = 1:size(allKeys,2)
                    hitKey = thiskey(istim);
                    thisstim = allstims(istim);
                    clusterinds = find(IDX == hitKey);
                    hitinds = find(IDXhatsub(clusterinds)==thisstim);
                    hits(ikey,istim,iestimate) = size(hitinds,1);
                    stimacc(ikey,istim,iestimate) = hits(ikey,istim,iestimate)/size(clusterinds,1);
                end
                totalkeyacc(ikey,iestimate) = sum(hits(ikey,:,iestimate))/size(IDX,1);
                avgkeyacc(ikey,iestimate) = mean(stimacc(ikey,:,iestimate));
            end
            
            
        end
        MaxTotKeyacc{iexpt,ib} = max(totalkeyacc);% max(mean(totalkeyacc,2));
        MaxAvgKeyacc{iexpt,ib} = max(avgkeyacc); %max(mean(avgkeyacc,2));
        %         MaxStimacc(ib) = mean(max(max(stimacc)));
    end
    %     end
%     matname = 'KmeansAccSeriesVmMeta5stims';
    matname = 'KmeansAccMeanVmMeta5stims';
%     matname = 'KmeansAccSeriesVmMeta5stimsBaseline';
%     matname = 'KmeansAccMeanVmMeta5stimsBaseline';
    if getnoisefloor == 1
        matname = 'KmeansAccVmMeta5StimsNoise';
    end
    save([r.Dir.Expt '/Analysis/KmeansAccVm/' matname '.mat'], ...
        'MaxTotKeyacc', 'MaxAvgKeyacc','binsizeVec')
end

%%
matname = 'KmeansAccMeanVmMeta5Stims';
load([r.Dir.Expt 'Analysis/KmeansAccVm/' matname '.mat']);

figure;
hold on
color = [0.5 0.5 0.5];
if getnoisefloor == 1;
    color = 'r';
end
for iexpt = 1:size(MaxTotKeyacc,1)
    subplot(5,5,iexpt)
    hold on
    matname = 'KmeansAccMeanVmMeta5Stims';
    load([r.Dir.Expt 'Analysis/KmeansAccVm/' matname '.mat']);
    for ib = 1:size(MaxTotKeyacc,2)
        if ~isempty(MaxTotKeyacc{iexpt,ib})
            color = 'b';
            scatter(binsizeVec(ib)*expt.wc.dt,mean(MaxTotKeyacc{iexpt,ib}),100,color,'fill')
        end
    end
    matname = 'KmeansAccMeanVmMeta5StimsBaseline';
    load([r.Dir.Expt 'Analysis/KmeansAccVm/' matname '.mat']);
    for ib = 1:size(MaxTotKeyacc,2)
        if ~isempty(MaxTotKeyacc{iexpt,ib})
            color = 'r';
            scatter(binsizeVec(ib)*expt.wc.dt,mean(MaxTotKeyacc{iexpt,ib}),100,color,'fill')
        end
    end
    ylabel('mean max acc','FontSize',14)
    xlabel(['binsize(s)  ' repexpts{iexpt}],'FontSize',14,'Interpreter','none')
    set(gca,'YLim',[0,1],'XLim',[0,1])
    
end
% ylabel('max (mean across subsamples) total hits rate','FontSize',14)
% xlabel('binsize(seconds)','FontSize',14)
% set(gca,'YLim',[0,1])
% title([num2str(numsubsamp) 'stimuli   ' num2str(numrep) 'estimates'], 'FontSize',14)


%%

figure;
set(gcf,'Color',[1,1,1])
hold on

matname = 'KmeansAccSeriesVmMeta5stims';
load([r.Dir.Expt 'Analysis/KmeansAccVm/' matname '.mat']);

color = 'm';%[0.5, 0.5, 0.5];%'k';%[0.8, 0.8, 0.8];%[0.3 0.3 0.3];
if getnoisefloor == 1;
    color = 'r';
end

binmag = [];
for ib = 1:size(MaxTotKeyacc,2)
    mag = [];
    exptind = 1;
    for iexpt = 4 %1:size(MaxTotKeyacc,1)
        if isempty(MaxTotKeyacc{iexpt,ib})
            mag(exptind) = nan;
            exptind = exptind +1;
        end
        if ~isempty(MaxTotKeyacc{iexpt,ib})
            mag(exptind) = mean(MaxTotKeyacc{iexpt,ib});
%             scatter(binsizeVec(ib)*expt.wc.dt,mean(MaxTotKeyacc{iexpt,ib}),100,color,'fill')
            exptind = exptind +1;
        end
    end
    binmag(ib) = nanmean(mag);
end
scatter(binsizeVec(1:size(MaxTotKeyacc,2))*expt.wc.dt,binmag,100,color,'fill')
ylabel('clustering hit rate','FontSize',18)
xlabel('window size (seconds)','FontSize',18)
set(gca,'YLim',[0,1])
title([num2str(numsubsamp) 'stimuli   ' num2str(numrep) 'estimates'], 'FontSize',14)

%
% matname = 'KmeansAccVm';
% if getnoisefloor == 1
%     matname = 'KmeansAccVmNoise';
% end
% save([r.Dir.Expt '/Analysis/KmeansAccVm/' matname '.mat'], ...
%     'MaxTotKeyacc', 'MaxAvgKeyacc','binsizeVec')
%%
matname = 'KmeansAccVmNoise';
load([r.Dir.Expt '/Analysis/KmeansAccVm/' matname '.mat']);

color = 'b';

figure;
hold on
matname = 'KmeansAccVm';
% load([r.Dir.Expt '/Analysis/KmeansAccVm/' matname '.mat']);

scatter(binsizeVec*expt.wc.dt,MaxTotKeyacc,30,'b','fill')

% matname = 'KmeansAccVmNoisetmp';
% load([r.Dir.Expt '/Analysis/KmeansAccVm/' matname '.mat']);
% for iestimate = 1:numestimates
% scatter(binsizeVec*expt.wc.dt,MaxTotKeyacc(:,iestimate),30,'r','fill')
% end
ylabel('max (mean across subsamples) total hits rate','FontSize',14)
xlabel('binsize(seconds)','FontSize',14)


figure;
hold on
matname = 'KmeansAccVm';
load([r.Dir.Expt '/Analysis/KmeansAccVm/' matname '.mat']);
for iestimate = 1:numestimates
    scatter(binsizeVec*expt.wc.dt,MaxAvgKeyacc(:,iestimate),30,'b','fill')
end
matname = 'KmeansAccVmNoisetmp';
load([r.Dir.Expt '/Analysis/KmeansAccVm/' matname '.mat']);
for iestimate = 1:numestimates
    scatter(binsizeVec*expt.wc.dt,MaxAvgKeyacc(:,iestimate),30,'r','fill')
end
ylabel('max (mean across subsamples) avg stimulus hits rate','Fontsize',14)
xlabel('binsize(seconds)','FontSize',14)


figure;
hold on
scatter(binsizeVec*expt.wc.dt,MaxStimacc,30,color,'fill')
ylabel('overall maximum (mean across subsamples) stimulus hit rate','FontSize',14)
xlabel('binsize(seconds)','FontSize',14)



%%
save([r.Dir.Expt '/Analysis/KmeansAccVm/KmeansAccVm.mat'],'MaxKeyacc','MaxStimacc','numBins','numUnique','numUniqueTemplate')

%%
exptnames = {
    'KP_B130_131120_p1c2a.mat'};
numbinind = 1;
numUnqAccEst = [];
sil_hat = [];

for ib = 1:size(binsizeVec,2)
    binsize = binsizeVec(ib);
    IDXhat = [];
    Chat = [];
    X = [];
    Xavg = [];
    xind = 1;
    for iexpt = 1:size(exptnames,1)
        thisexpt = exptnames{iexpt};
        load([r.Dir.Expt thisexpt]); %load the experiment .mat file
        vmexpt = filtesweeps(expt,0,'Vm',0); %filter out all trials that were not recorded under current clamp at 0 holding
        stimcond = expt.stimcond; %the structure with information about the stimuli played
        for istim = 1:size(stimcond,2) %grab only loudest db of db stims... or non-db stims from this set of stimuli
            notdbstim = [];
            if isempty(regexp(stimcond(istim).wavnames,'d'))
                notdbstim = [notdbstim, istim];
            end
            loudstim = [];
            if ~isempty(regexp(stimcond(istim).wavnames,'d'))
                if ~isempty(regexp(stimcond(istim).wavnames,'d80'))
                    loudstim = [loudstim,istim];
                end
            end
        end
        
        stiminds = union(notdbstim, loudstim); %left with only "unique" stimuli....
        
        if isempty(stiminds)
            continue
        end
        %bin the response for that stimulus into substimuli of "binsize" long
        for istim = 1:max(size(stiminds))
            thisind = stiminds(istim);
            stimwav = stimcond(thisind).wavs; %this is the stimulus wavform in case you want to see it
            %(yes... it is a repeat of a motif, but that could end up
            %being something interesting to note... because the
            %response during the second repeat is usually different
            %from the first repeat even though the stimulus is exactly
            %the same... context/history effects
            stimexpt = filtesweeps(vmexpt,0,'wavnames',stimcond(thisind).wavnames); %filtering out all trials that is not that stimulus
            sigon = expt.analysis.params.baselinewin(2);
            sigoff = sigon+round(size(stimcond(thisind).wavs,1)/44100/expt.wc.dt);
            sigdata = medfilt1(stimexpt.wc.data(:,sigon:sigoff),[],200,2); %all trials in response to that stimulus
            ntrials = size(sigdata,1);
            numbins = floor((sigoff-sigon)/binsize); %how many subresponses of "binsize" width can there be
            bins = [1:binsize:numbins*binsize]; %make the edges of the subresponse bins
            for ibin = 1:numbins-1 %parse the response into subresponses
                X = [X;sigdata(:,bins(ibin):(bins(ibin+1)-1))]; %stack all of the trials of all of the subresponses on top of each other
                IDXhat = [IDXhat;repmat(xind,size(sigdata,1),1)]; %keep track of what bin each trial "belongs to"
                [idx,c] = kmeans(sigdata(:,bins(ibin):(bins(ibin+1)-1)),1); %get an estimate of the "real" cluster center
                Chat (ibin,:) = c; %keep track of the centroid for each set of trials for each bin
                Xavg = [Xavg;mean(sigdata(:,bins(ibin):(bins(ibin+1)-1)))]; %this is the mean Vm response for each bin
                xind = xind +1;
            end
        end
    end
    
    %I added this "template" part in as an idea for maybe how to constrain the
    %clustering a little so that it was more specific, but I am not sure that
    %it helps and ultimately it might be too non-agnostic to the trials data (X matrix clustering) to use it
    
    
    %     ik = numbins-1;
    ik = numUnique(ib);
    [IDXtemplate, Ctemplate] = kmeans(Xavg,ik,'replicate',5); %calculate the centroid of the average Vm response for each subresponse
    %     [ix,b] = sort(IDXtemplate)
    %     IDXtemplate = IDXtemplate(b);
    %     Ctemplate = Ctemplate(b,:);
    %
    
    
    replicates = 10; %how many times to run kmeans
    
    [IDX, C, sumd,D] = kmeans(X,ik,'replicate',replicates); %,'start',...
    %repmat(Ctemplate,[1,1,replicates]));
    
    [silh3,h] = silhouette(X,IDX);
    sil_hat{ib} = silh3;
    %results of kmeans... IDX is the indices of what cluster it belongs to... and C is the centroid locations
    %Ctemplate was used to set the starting points for the centroids
    %for the analsysis...
    
    
    
    IDXstims = reshape(IDX,ntrials,size(IDX,1)/ntrials);
    maxstim = numBins(ib);
    repeatedbins = maxstim - ik;
    scatteredtrials = 0;
    if repeatedbins ~= 0
        scatteredtrials = ceil(repeatedbins/ik);
    end
    accthresh(ib) = 6/(ntrials+scatteredtrials);
    
    maxID = [];
    acc = [];
    for istim = 1:size(IDXstims,2)
        thisstim = IDXstims(:,istim);
        unqclus = unique(thisstim);
        numclus = [];
        for iclus = 1:size(unqclus,1)
            numclus(iclus) = size(find(thisstim == unqclus(iclus)),1);
        end
        nummax = size(find(numclus == max(numclus)),2);
        if nummax == 1
            maxID(istim) = unqclus(find(numclus == max(numclus)));
            acc(istim) = size(find(thisstim == maxID(istim)),1)/size(IDXstims,1);
            
            numstim = ik; % = numUnique(ib);
            %what is the "best accuracy" that can be achieved given the
            %#clusters?
            
        end
        if nummax >= 2
            maxID(istim) = nan;
            acc(istim) = nan;
        end
    end
    
    isunq = [];
    unqID = unique(maxID(~isnan(maxID)));
    for iunq = 1:size(unqID,2)
        if size(find(maxID == unqID(iunq)),2) > 1
            isunq(iunq) = 0;
        end
        if size(find(maxID == unqID(iunq)),2) == 1
            useClus = find(maxID == unqID(iunq));
            if acc(useClus) < accthresh(ib)
                isunq(iunq) = 0;
            end
            if acc(useClus) >= accthresh(ib)
                isunq(iunq) = 1;
            end
        end
    end
    numUnqAccEst(ib) = size(find(isunq),2);
    
    
    
    
    %get distance from each estimated cluster to every "real cluster"
    %     D = zeros(size(C,1),size(Ctemplate,1));
    %     nclusts = size(C,1);
    %     %'sqeuclidean'
    %     for i = 1:size(C,1)
    %         D(:,i) = (C(:,1) - Ctemplate(i,1)).^2;
    %         for j = 2:size(C,2)
    %             D(:,i) = D(:,i) + (C(:,j) - Ctemplate(i,j)).^2;
    %         end
    %     end
    %
    %     %need to figure out more how to match these up...
    %     minD_Ctemp = [];
    %     IDXpredict = [];
    %     for iclus = 1:size(D,2)
    %         minD_Ctemp(iclus) = find(D(:,iclus) == min(D(:,iclus)));
    %     end
    %     for iclus = 1:size(D,2)
    %         IDXpredict = [IDXpredict; repmat(minD_Ctemp(IDXtemplate(iclus)),ntrials,1)];
    %     end
    %
    %     prediction = [];
    %     for iresp = 1:size(IDXpredict,1)
    %         if IDXpredict(iresp) == IDX(iresp)
    %             prediction(iresp) = 1;
    %         end
    %         if IDXpredict(iresp) ~= IDX(iresp)
    %             prediction(iresp) = 0;
    %         end
    %     end
    %
    %    predictioncorrect(ib) = size(find(prediction),2)/size(prediction,2);
    %
    % %      figure;hold on
    % %      sigdata_template = [];
    % %      Ctempinds = [1:numbins-1];
    % %     for iclus = 1:size(Ctemplate,1)
    % %         sigdata_template{iclus} = Xavg(find(IDXtemplate == IDXtemplate(iclus)),:);
    % %         subplot(size(Ctemplate,1),1,iclus)
    % %         plot((sigdata_template{iclus}));
    % %     end
    %
    %     %this should be the same as above
    %     if doplot ==1
    %         figure;hold on
    %         sigdata_resp = [];
    %         for iclus = 1:size(Ctemplate,1)
    %             sigdata_resp{iclus} = X(find(IDXhat == iclus),:);
    %             subplot(size(Ctemplate,1),1,iclus)
    %             plot(mean(sigdata_resp{iclus}));
    %         end
    %
    %         figure;hold on
    %         sigdata_hat = [];
    %         Cinds = [1:numbins-1];
    %         for iclus = 1:size(D,2)
    %             sigdata_hat{iclus} = X(find(IDX == (minD_Ctemp(IDXtemplate(iclus)))),:);
    %             subplot(size(D,2),1,iclus)
    %             plot(mean(sigdata_hat{iclus}));
    %             %title each as its corresponding cluster ID
    %             title(['cluster ' num2str(minD_Ctemp(IDXtemplate(iclus))) ...
    %                 ' similar to response # '  num2str(iclus) ' from template '])
    %         end
    %     end
    
end
figure;scatter(binsizeVec*expt.wc.dt,numUnqAccEst./numUnique)
hold on
line(binsizeVec*expt.wc.dt,accthresh)

figure;
hold on
line(binsizeVec*expt.wc.dt,numBins,'color','g')
line(binsizeVec*expt.wc.dt,numUnqAccEst,'color','k')
legend({'# bins','# accurate'})
%%
numUnique(ib)

bintime = binsizeVec*expt.wc.dt;
figure;
scatter(bintime,predictioncorrect,50,'k','fill')
xlabel('binsize in milliseconds','FontSize',14)
ylabel('percent "accurate" classification by kmeans','FontSize',14)

empsil_x=[];
empsil_p=[];
erbarsil = [];
medsil = [];
for isil = 1:size(sil_hat,2)
    thissil = sil_hat{isil};
    medsil(isil) = median(thissil);
    Q  =  quantile(thissil,[0.15,0.85]);
    erbarsil(isil,:) = Q;
    [x,p ]= empcdf(thissil);
    empsil_x{isil} = x';
    empsil_p{isil} = p';
end


figure
hold on
scatter(bintime, medsil)
for isil = 1:size(sil_hat,2)
    errorbar(bintime(isil),medsil(isil),(medsil(isil)-erbarsil(isil,1)), ...
        (erbarsil(isil,2)-medsil(isil)))
end
xlabel('binsize in milliseconds','FontSize',14)
ylabel('median sillouhette values +- 85%','FontSize',14)

figure
hold on
[grad,im]=colorGradient([0.8, 0.8, 0.8],[0, 0, 0],size(sil_hat,2));
for isil = 1:size(sil_hat,2)
    stairs(empsil_x{isil},empsil_p{isil},'color',grad(isil,:))
    
end
xlabel('sillouhette values across all trials','FontSize',14)
ylabel('cdf for each binsize (narrow:gray to wide:black)','FontSize',14)

%%
figure %plot the sillouettes of the clusters to see how "good" the clustering is
[silh3,h] = silhouette(X,IDX);
set(get(gca,'Children'),'FaceColor',[.8 .8 1])
xlabel('Silhouette Value')
ylabel('Cluster')


for isort = 1:numbins-1
    numclus(isort) = size(find(IDX == isort),1);
end
numclus %seems like I want to figure a way to iterate this process so that the clusters are even
%sized... even if it makes the clusters "worse." Because the question is
%whether the trials were appropriately assigned given certain constraints
%(like there were 10 trials per stimulus and there are numbins-1 stimuli)

IDXrevise = IDX;
for iclus = 1:size(numclus,2)
    maxind = find(numclus == max(numclus) )
    clusinds = find(IDX == maxind);
    silh3_clus = silh3(clusinds);
    [b,ix] = sort(silh3_clus);
    notclusinds = clusinds(ix(11:end));
    allclus = [1:numbins-1];
    otherclus = allclus(find(allclus~=maxind));
    for inot = 1:size(notclusinds,2) %for each point that is furthest from centroid... find nearest neighbor
        otherD = D(notclusinds(inot),:);
        [b,ix] = sort(otherD);
        switchclus = ix(2);
        IDXrevise(notclusinds(inot)) = switchclus;
    end
    numclus = [];
    for isort = 1:numbins-1
        numclus(isort) = size(find(IDXrevise == isort),1);
    end
    numclus
    
    for ii = 1:size(numclus,2) %get new centroid matrix
        
    end
end


%DISTFUN Calculate point to cluster centroid distances.
[n,p] = size(X);
D = zeros(n,size(C,1));
nclusts = size(C,1);
%'sqeuclidean'
for i = 1:nclusts
    D(:,i) = (X(:,1) - C(i,1)).^2;
    for j = 2:p
        D(:,i) = D(:,i) + (X(:,j) - C(i,j)).^2;
    end
    % D(:,i) = sum((X - C(repmat(i,n,1),:)).^2, 2);
end


%------------------------------------------------------------------

%function [centroids, counts] = gcentroids(X, index, clusts, dist)
%GCENTROIDS Centroids and counts stratified by group.
p = size(X,2);
num = length(clusts);
centroids = NaN(num,p);
counts = zeros(num,1);

for i = 1:num
    members = (index == clusts(i));
    if any(members)
        counts(i) = sum(members);
        switch dist
            case 'sqeuclidean'
                centroids(i,:) = sum(X(members,:),1) / counts(i);
            case 'cityblock'
                % Separate out sorted coords for points in i'th cluster,
                % and use to compute a fast median, component-wise
                Xsorted = sort(X(members,:),1);
                nn = floor(.5*counts(i));
                if mod(counts(i),2) == 0
                    centroids(i,:) = .5 * (Xsorted(nn,:) + Xsorted(nn+1,:));
                else
                    centroids(i,:) = Xsorted(nn+1,:);
                end
            case {'cosine','correlation'}
                centroids(i,:) = sum(X(members,:),1) / counts(i); % unnormalized
            case 'hamming'
                % Compute a fast median for binary data, component-wise
                centroids(i,:) = .5*sign(2*sum(X(members,:), 1) - counts(i)) + .5;
        end
    end
end
% function