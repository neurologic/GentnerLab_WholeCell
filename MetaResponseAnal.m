trials = 5;
signalDB = 80;
nresid_stim = [];
nresid_base = [];
s_base = [];
s_stim = [];
stimind = 1;
for iexpt=1:size(repexpts,2)
    iexpt
    thisexpt=repexpts{iexpt};
    load([r.Dir.Expt thisexpt])
    
    vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
    
    table=getClampTab(expt,{'clamp',0});
    highpassdata=HighpassGeneral(vmexpt.wc.data,[],100,1/expt.wc.dt);
    
    %     outcell=PlotAndAsk(highpassdata,'spikesthresh','negative');
    %     cellfun(@eval,outcell);
    % for these "10rep expts", i know that all the spikesthresh are 0.01 at
    % this point... so hard code for now
    spikesthresh = 0.01;
    negative = 0;
    
    keepsigs=reprequire(table,trials);
    allsig=table.sigsplayed;
    stimcond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
    repsigs=allsig(keepsigs);
    basetimes=basewin{iexpt};
    
    [sigon,sigoff]=GetSigTimes(expt,expt.stimcond,1);
    
    [dbstimcond,dblevels]=getDBstimcond(vmexpt);
    if isempty(dbstimcond)
        continue
    end
    
    for istimcond=1:size(dbstimcond,2)
        thiscond=dbstimcond{istimcond};
        thisdb=dblevels{istimcond};
        useind=[];
        for idb=1:size(thisdb,2)
            thiscond(idb).wavnames;
            testrep=regexp(repsigs,thiscond(idb).wavnames);
            isrep=0;
            for itest=1:size(testrep,1)
                if ~isempty(testrep{itest})
                    isrep=1;
                end
            end
            if isrep==1
                useind(idb)=1;
            else useind(idb)=0;
            end
        end
        if size(find(useind),2)<2
            continue
        end
        thiscond=thiscond(find(useind));
        thisdb=thisdb(find(useind));
        
        %         for now just use the 80dB stim
        dbind = find(thisdb==signalDB);
        if isempty(dbind)
            continue
        end
        sigexpt = filtesweeps(vmexpt,0,'wavnames',...
            thiscond(dbind).wavnames);
        sigdata = sigexpt.wc.data*1000;
        sigdata_filt = medfilt1(sigdata,200,[],2);
        
        %get residuals
        resid_data = sigdata_filt - repmat(mean(sigdata_filt),size(sigdata_filt,1),1);
        resid_edges = [-40:2:40];
        tmp_stim = resid_data(:,sigon:sigoff);
        tmp_stim = reshape(tmp_stim',1,...
            size(tmp_stim,1)*size(tmp_stim,2));
        tmp_base = resid_data(:,basetimes(1):sigon);
        tmp_base = reshape(tmp_base',1,...
            size(tmp_base,1)*size(tmp_base,2));
        nresid_stim(stimind,:) = histc(tmp_stim,resid_edges)./size(tmp_stim,2);
        nresid_base(stimind,:) = histc(tmp_base,resid_edges)./size(tmp_base,2);
        
        s_base (stimind) = skewness(tmp_base);
        s_stim (stimind) = skewness(tmp_stim);
        
%         figure;
%         hold on
%         stairs(resid_edges,nresid_base(stimind,:),'color','b','LineWidth',3)
%         stairs(resid_edges,nresid_stim(stimind,:),'color','r','LineWidth',3)
%         scatter(mean(tmp_stim),max([nresid_stim(stimind,:),nresid_base(stimind,:)])+0.05,100,'r','v')
%         scatter(mean(tmp_base),max([nresid_stim(stimind,:),nresid_base(stimind,:)])+0.05,100,'b','v')
%         ylims = get(gca,'YLim');
%         text(resid_edges(10),ylims(2)-0.05,{['base skew = ' num2str(s_base(stimind))];...
%             ['stim skew = ' num2str(s_stim(stimind))]},'HorizontalAlignment','center');
%         %using residuals is a really easy way to "highpass"
%         title([expt.name ';  stimulus# ' num2str(istimcond)...
%             ';  distribution Vm residuals;  stim = ' ...
%             num2str(thisdb(dbind)) 'dBSPL'],'Interpreter','none');
%         
        stimind = stimind +1;
        
        %probability of a spike for a give Vm
        %during baseline vs during stimulus
        %subtract residuals and get highpass data for 
%%%%% based on how some of these residuals look, need to eliminate
%%%%% some trials and/or some cell:signal pairs?
%         highpassdata = sigdata - repmat(mean(sigdata_filt),size(sigdata,1),1);
%         figure;plot(highpassdata')
%         threshold = input('what is spike threshold?');
         highpassdata=HighpassGeneral(sigdata,[],metadata.highcutoff,1/expt.wc.dt); 
         threshold = 20;
        [spikesmat, gausstosmooth]=getspikesmat(highpassdata,threshold,expt);
        spiketimes = [];
        for itrial=1:size(spikesmat,1)
            spiketimes{itrial}=find(spikesmat(itrial,:));
        end
        
        vm_edges = [-80:2:-30];
        x_bins = vm_edges(1,1:end-1) + mean(diff(vm_edges))/2;
        p_spk_base = [];
        p_spk_resp = [];
        for ibin = 1:size(vm_edges,2) - 1;
            thisbin = [vm_edges(ibin),vm_edges(ibin+1)];
            numsamps_resp = [];
            numsamps_base = [];
            spks_resp = [];
            spks_base = [];
            
            for itrial=1:size(spikesmat,1)
                spks_trial = spiketimes{itrial};
                % p(spike:Vm) during the baseline on trial by trial
                resp_trial = sigdata_filt(itrial,sigon:sigoff);
                sect1 = find(resp_trial >= thisbin(1));
                sect2 = find(resp_trial < thisbin(2));
                intersectinds = intersect(sect1,sect2)+sigon; %need to make these inds match spiketimes
                if ~isempty(intersectinds)
                    numsamps_resp = size(numsamps_resp,2) + size(intersectinds,2);
                    spks_resp = size(spks_resp,2) + size(intersect(intersectinds,spks_trial),2);
                end
                % p(spike:Vm) during stimulus on trial by trial
                base_trial = sigdata_filt(itrial,basetimes(1):sigon);
                sect1 = find(base_trial >= thisbin(1));
                sect2 = find(base_trial < thisbin(2));
                intersectinds = intersect(sect1,sect2)+basetimes(1);%need to make these inds match spiketimes
                if ~isempty(intersectinds)
                    numsamps_base = size(spks_base,2) + size(intersectinds,2);
                    spks_base = size(spks_base,2) + size(intersect(intersectinds,spks_trial),2);
                end
            end
            
            if isempty(numsamps_resp)
                p_spk_resp (ibin) = NaN;
            elseif ~isempty(numsamps_resp)
                p_spk_resp (ibin) = spks_resp/numsamps_resp;
            end
            
            if isempty(numsamps_base)
                p_spk_base(ibin) = NaN;
            elseif ~isempty(numsamps_base)
                p_spk_base (ibin) = spks_base/numsamps_base;
            end
        end
        max_pspk = max([max(p_spk_base),max(p_spk_resp)]);
        P_spk_base(stimind,:) = p_spk_base / max_pspk;
        P_spk_resp(stimind,:) = p_spk_resp / max_pspk;
        
        stimind = stimind + 1;
    end
    
end

figure;
hold on
stairs(resid_edges,mean(nresid_base),'color','b','LineWidth',3)
stairs(resid_edges,mean(nresid_stim),'color','r','LineWidth',3)
axis tight
ylims = get(gca,'YLim');
text(resid_edges(10),ylims(2)-0.05,{['mean base skew = ' num2str(mean(s_base))];...
    ['mean stim skew = ' num2str(mean(s_stim))]},'HorizontalAlignment','center');
title([num2str(stimind-1) ' cell:signal pairs    stimulus at ' ...
    num2str(signalDB) 'dBSPL'])

%using residuals is a really easy way to "highpass"

%p(spk|Vm)
figure;
hold on
scatter(x_bins,nanmean(P_spk_base,1),100,'b','fill')
scatter(x_bins,nanmean(P_spk_resp,1),100,'r','fill')
xlabel('Vm')
ylabel('normalized p(spike|Vm)')
set(gca,'XTick', x_bins, 'XTickLabel', x_bins)
title([num2str(stimind-1) ' cell:signal pairs    stimulus at ' ...
    num2str(signalDB) 'dBSPL'])
