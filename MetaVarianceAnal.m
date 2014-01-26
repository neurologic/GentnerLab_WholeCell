%%
r=rigdef('mac')
%%


save([r.Dir.Expt 'Analysis/Meta_10Trial_Responses_DB/exptsForICanal.mat'],...
    'exptnames','readyexpt','sigdata','vardata','basewin','allstims','repexpts',...
    'addstimcond','vmresp','base_vmresp','base_spkresp','spkresp','respDB',...
    'respdata','cb','cr','respslope','baseslope')

%%
load([r.Dir.Expt 'Analysis/Meta_10Trial_Responses_DB/exptsForICanal.mat'],...
    'exptnames','readyexpt','sigdata','vardata','basewin','allstims','repexpts',...
    'addstimcond','vmresp','base_vmresp','base_spkresp','spkresp','respDB',...
    'respdata','cb','cr','respslope','baseslope')
%%
d=dir(r.Dir.Expt);
exptnames=[];
exptind=1;
for id=1:length(d)
    thisname=d(id).name;
    if ~isempty(regexp(d(id).name,'.mat'))
        if ~isempty(regexp(d(id).name,'KP'))
            exptnames{exptind}=d(id).name;
            exptind=exptind+1;
        end
    end
end
addstimcond=[];
addind=1;
%%
sigdata = [];
vardata = [];
basewin = [];
allstims = [];
repexpts = [];
trials=8;
exptind=1;
%%
for iexpt=1:size(readyexpt,1)
    
    thisexpt=readyexpt{iexpt};
    %load in expt
    load([r.Dir.Expt thisexpt])
    
    %check that the expt has a folder
    isfold=exist([r.Dir.Expt expt.name], 'dir');
    if isfold~=7
        %make the folder if it didn't exist
        mkdir(r.Dir.Expt,expt.name);
    end
    
    %check if this expt has a stimcond field
    hasstim=isfield(expt,'stimcond');
    if ~hasstim
        %enter this expt in a list to go through later and add stimcond
        %skips this expt for now
        addstimcond{addind}=expt.name;
        addind=addind+1;
        continue
    end
    
    %find if this expt has a table
    hastable=isfield(expt,'table');
    if ~hastable
        %query experiment if it does not have a table yet
        [expt,table,hfig]=QueryExpt(expt);
    end
    
    %find if this expt had IC trials
    %find if this expt had IC trials
    [isclamp,icID]=findclamp(expt,'ic');
    if isclamp==1
        %if there is more than one IC trial, pick the 0 holding current
        if size(icID,2)>1
            useVm=[];
            clamps=[];
            for iclamp=1:size(icID,2)
                clamps(iclamp)=expt.table(icID(iclamp)).clamp;
            end
            if ~isempty(find(clamps==0))
                useVm=0;
            end
            %if no zero holding current, use the next most negative
            if isempty(find(clamps==0))
                clamps=sort(clamps);
                [s,t]=crossing(clamps);
                useVm=clamps(s);
            end
        end
        if size(icID,2)==1
            %if only one holding current use that one
            useVm=expt.table(icID).clamp;
        end
        
        %mark this expt as "used" for this analysis
        
        table=getClampTab(expt,{'clamp',useVm});
        keepsigs=reprequire(table,trials);
        stimcond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
        
        if ~isempty(stimcond)
            usedexpt(iexpt)=1;
            
            vmexpt=filtesweeps(expt,0,'Vm',useVm);
            
            %check the waveonset_time and baselinewin
            %add .check to expt.analysis.params if it is checked
            checked_params = isfield(expt.analysis.params,'checked');
            if ~checked_params
                basetimes=vmexpt.analysis.params.baselinewin;
                s=PlotAndAsk(vmexpt.wc.data(:,basetimes(1):basetimes(2)),'basebegin')
                cellfun(@eval,s);
                if ~isempty(basebegin)
                    expt.params.baselinewin(1)=basetimes(1)+basebegin;
                end
                hfig = figure;
                hold on
                hl = line([1:size(expt.wc.data,2)]*expt.wc.dt,vmexpt.wc.data');
                [sigon,sigoff]=GetSigTimes(vmexpt,stimcond,1);
                
                SigTimeBox(gca, sigon*expt.wc.dt, sigoff*expt.wc.dt, get(gca,'YLim'),'k')
                onset_correct = input('sigonset correct?');  % 1:true, 0:false
                if ~onset_correct
                    onset = input('what is it instead?'); % in second...
                    %(because will be either 1.3 or 1.8... haven't done any
                    %other times in my time here)
                    expt.analysis.params.waveonset_time = onset/expt.wc.dt;
                end
                expt.analysis.params.checked = [1];
                save(fullfile(r.Dir.Expt,expt.name),'expt')
                
            end
            
            [s,v] = GetVarDistrib(vmexpt,stimcond);
            sigdata{exptind} = s;
            vardata{exptind} = v;
            allstims{exptind} = stimcond;
            repexpts{exptind} = thisexpt;
            basewin{exptind} = expt.analysis.params.baselinewin;
            exptind=exptind+1;
            
        end
        if isempty(stimcond)
            usedexpt(iexpt)=0;
            continue
        end

    end
    if isclamp==0
        usedexpt(iexpt)=0;
    end
    
    %mark this expt as "used" for this analysis
    
    
    
end
%% did this to get ready expt... but now readyexpt is saved.
%all the rest of the expts need to have stimcond added, etc to be able to
%be filtered and used with this workflow
reind=find(usedexpt);
readyexpt=exptnames(reind)

%% getting Rs and Rin to find 
for iexpt=1:length(readyexpt)
    thisexpt=readyexpt{iexpt};
    %load in expt
    load([r.Dir.Expt thisexpt])
    vmexpt=filtesweeps(expt,0,'Vm',0);
    outcell=PlotAndAsk(stepdata,['stepon','stepoff']);
    xtime=[553:1000];
    stepdata=nanmean(vmexpt.wc.data(:,553:1000))-max(nanmean(vmexpt.wc.data(:,553:1000)));
    f=fit([553:1000]',stepdata','exp2')
    %     f=fit(xtime',nanmeanstepdata','exp2');
    figure; hold on
    
    a=f.a;
    b=f.b;
    c=f.c;
    d=f.d;
    alltau= [f.b,f.d];
    alltau=alltau(find(alltau<0));
    tau=round(nanmean(alltau));
    xtimestep=[1:size(stepdata,2)]*expt.wc.dt;
    line(xtimestep,nanmean(stepdata,1),'color','k','LineWidth',2)
    line(xtimestep,(c*exp(tau*xtimestep))+maxstepdata+a);
    minVe=min(min(xtimestep,(c*exp(d*xtimestep))+maxstepdata+a));
    maxVe=max(max(nanmean(stepdataNorm,1)));
    deltaVe=maxVe-minVe;
    stepamp=-0.075;
    Rs=deltaVe/0.075;
    
    rs(iexpt)=expt.table.Rs
end

%%
% 
% %take note of which stimuli are DB stimuli so can sort out that data
% isdb=[];
% thisdb=[];
% for istimcond=1:size(allstims,2)
%     thiscond=allstims{istimcond};
%     for istim=1:size(thiscond,2)
%         thisstim=thiscond(istim).wavnames;
%         dind=regexp(thisstim,'d');
%         if ~isempty(dind)
%             isdb{istimcond}(istim)=1;
%             thisdb{istimcond}(istim)=str2num(thisstim(1,end-1:end));
%         end
%         if isempty(dind)
%             isdb{istimcond}(istim)=0;
%         end
%     end
% end
% 
% dbsigdata=[];
% dbvardata=[];
% dbon=[];
% dboff=[];
% dbpercsig=[];
% percbase=[];
% dballstims=[];
% 
% dbins=[25:20:85];
% %for each of the cells, separate data for each of the db categories
% for idb=1:size(dbins,2)-1
%     needcross=idb;
%     stimind=1;
%     for istimcond=1:size(allstims,2)
%         celldb=thisdb{istimcond};
%         for istim=1:size(celldb,2)
%             if crossing(dbins,[],celldb(istim))==needcross
%                 dbon (idb,stimind)= onoff{istimcond}(istim,1);
%                 dboff(idb,stimind) = onoff{istimcond}(istim,2);
%                 dbpercbase(idb,stimind)=percsignif{istimcond}(istim,1);
%                 dbpercsig(idb,stimind)=percsignif{istimcond}(istim,2);
%                 dbvardata{idb,stimind}=vardata{istimcond}(istim,:);
%                 dbsigdata{idb,stimind}=sigdata{istimcond}(istim,:);
%                 stimind=stimind+1;
%             end
%         end
%     end
%     oldind(idb)=stimind;
% end


%%
isdb=[];
thisdb=[];
allsigdata=[];
allvardata=[];
allbasewin=[];
stimind=1;
for istimcond=1:size(allstims,2)
    thiscond=allstims{istimcond};
    [sigon,sigoff]=GetSigTimes(expt,thiscond,1);
    
    for istim=1:size(thiscond,2)
        thisstim=thiscond(istim).wavnames;
        dind=regexp(thisstim,'d');
        if ~isempty(dind)
            isdb(stimind)=1;
            thisdb(stimind)=str2num(thisstim(1,end-1:end));
        end
        if isempty(dind)
            isdb(stimind)=NaN;
        end
        allsigdata{stimind}=sigdata{istimcond}(istim,:);
        allvardata{stimind}=vardata{istimcond}(istim,:);
        allbasewin{stimind}=basewin{istimcond};
        allsigwin{stimind}=[sigon,sigoff];
        allexptnames{stimind}=repexpts{istimcond};
        allsignames{stimind}=thisstim;
        stimind=stimind+1;
    end
end


%% this section is to get the variance in windows... 
% but it uses the arcane routine to get response windows from sigdata
sigind=1;
percresp=[];
percnot=[];
hbase2resp=[];
hbase2not=[];
medbound=[];
meanresp=[];
meannot=[];
avgVmResp=[];
avgVmNoresp=[];

% for iexpt=1:size(repexpts,1)-1;
%     base=basewin{iexpt};
%     for isig=1:size(sigdata{iexpt},1)
%         [xb,pb]=empcdf(sigdata{iexpt}(isig,base(1):sigon));
for isig=1:size(allsigdata,2);
    base_vmresp=allbasewin{isig};
    %find regions of significant Vm depol
    sigon=1.8/expt.wc.dt;
    sigoff=allsigwin{isig}(2);
    tmpdata=allsigdata{isig};
    [xb,pb]=empcdf(tmpdata(1,base_vmresp(1):sigon));
    tmp=find(pb>0.85);
    confbound=xb(min(tmp));
    respinds=find(tmpdata(sigon:sigoff)>confbound);
%     notrespinds=find(tmpdata(sigon:sigoff)<confbound);
    
    %find response windows ...THEN convert to respinds and notrespinds
    %that way I can window around the same responses I calculate changes in
    %dB for later... and I can make sure each response is actually
    %different and is long enough to count
    tmp=sort(respinds);
    t=unique(tmp); %across all db, all indices of sdata > confidence
    d=diff(t);
    %time between responses must be more than 50ms to count as
    %different resposnes
    difftime=round(0.0500/expt.wc.dt);
    respbreakind=find(d>difftime);
    respwin=min(t);
    for ibrk=1:size(respbreakind,2)
        respwin=[respwin,t(respbreakind(ibrk)),t(respbreakind(ibrk)+1)];
    end
    respwin=[respwin,max(t)];
    respwin=reshape(respwin',2,size(respwin,2)/2);
    respwin=respwin+sigon;
    resptimes=respwin(2,:)-respwin(1,:);
    %each response must last longer than 50ms to count as a response
    %window
    useresp=find(resptimes>round(0.0500/expt.wc.dt));
    respwin=respwin(:,useresp);
    %now convert to respinds full vector instead of window bounds
    respinds=[];
    for iresp=1:size(respwin,2)
        addinds=[respwin(1,iresp):respwin(2,iresp)];
        respinds=[respinds, addinds];
    end
    %now get notrespinds
    tmpvec=[1:size(tmpdata,2)];
    tmpvec(respinds)=0;
    notrespinds = find(tmpvec(1,sigon:sigoff)~=0)+sigon;
    
    %find regions of significantly low variance
    tmpvar=allvardata{isig};
    varallowt=round(0.1/expt.wc.dt); %use this filter to allow variance in real time to go above signif a little
    tmpvar=medfilt1(tmpvar,varallowt,[],2);
    ci = getCDFconf(tmpvar(1,base_vmresp(1):sigon),85);
    varinds=find(tmpvar(sigon:sigoff)<=ci(1));
    
    %test signif of distribution of var for resp and notresp periods
    %against baseline
    basevar=tmpvar(1,base_vmresp(1):sigon);
    
    if ~isempty(respinds)
        hbase2resp(sigind)=kstest2(basevar,tmpvar(respinds));
        %get percent of resp and notresp periods that has signif var
        numvar=intersect(varinds,respinds);
        percresp(sigind)=size(numvar,2)/size(respinds,2);
        %get "nanmean" variance during each period: baseline, response, not
        %response
        ntmpvar=tmpvar./max(tmpvar);
        meanresp(sigind)=nanmean(ntmpvar(respinds));
        avgVmResp(sigind)=mean(allsigdata{isig}(1,respinds));
    else
        hbase2resp  (sigind)=NaN;
        percresp(sigind)=NaN;
        meanresp(sigind)=NaN;
        avgVmResp(sigind)=NaN;
    end
    if ~isempty(notrespinds)
        hbase2not(sigind)=kstest2(basevar,tmpvar(notrespinds));
        numnot=intersect(varinds,notrespinds);
        percnot(sigind)=size(numnot,2)/size(notrespinds,2);
        ntmpvar=tmpvar./max(tmpvar);
        meannot(sigind)=nanmean(ntmpvar(notrespinds));
        avgVmNoresp(sigind)=mean(allsigdata{isig}(1,notrespinds));
    else
        hbase2not(sigind)=NaN;
        percnot(sigind)=NaN;
        meannot(sigind)=NaN;
        avgVmNoresp(sigind)=NaN;
    end
   
    
    
%     [xb,pb]=empcdf(ntmpvar(1,base(1):sigon));
%     medind=find(pb<=0.5);
%     medbound(sigind)=xb(max(medind));
    medbound(sigind)=nanmean(ntmpvar(1,base_vmresp(1):sigon));
   
    
    sigind=sigind+1;
    
end
%%
figure;
hold on
line([1:size(hbase2resp,2)],hbase2resp,'color','r','LineWidth',4);
line([1:size(hbase2not,2)],hbase2not,'color','b','LineWidth',2);

figure;
hold all
vecnm={'medbound','meanresp','meannot'};
for inm=1:size(vecnm,2)
    s=['[xb,pb]=empcdf(' vecnm{inm} ');'];
    eval(s);
    tmp=find(pb<0.85);
    upbound(inm)=xb(max(tmp));
    tmp=find(pb<0.15);
    lowbound(inm)=xb(max(tmp));
    stairs(xb,pb,'LineWidth',3)
end
legend(vecnm);

figure;
plotvec=[nanmean(medbound),nanmean(meanresp),nanmean(meannot)];
errorbar([0,1,2],plotvec,plotvec-lowbound,upbound-plotvec,'ok');
if kstest2(medbound,meanresp)
        line([0,1],[1,1],'color','r');
    end
    if kstest2(medbound,meannot)
        line([0,2],[1.2,1.2],'color','r');
    end
    if kstest2(meannot,meanresp)
        line([1,2],[1.4,1.4],'color','r');
    end
    set(gca,'YLim',[0,1.5])
     set(gca,'XTick',[0,1,2])
set(gca,'XTickLabel',{'spont','response','noresp'})
ylabel('average normalized variance ; ci=85%')
    
figure;
hold on
scatter(repmat(0,size(medbound,1),size(medbound,2)),medbound);
scatter(repmat(1,size(meanresp,1),size(meanresp,2)),meanresp);
scatter(repmat(2,size(meannot,1),size(meannot,2)),meannot);
  set(gca,'YLim',[0,1.5])
     set(gca,'XTick',[0,1,2])
set(gca,'XTickLabel',{'spont','response','noresp'})
ylabel('average normalized variance')
    
% need to figure out why the percent responses are not registering?
% so many of the entries are just zeroes
%
%
% figure;
% hold on
% scatter(repmat(1,size(percresp,1),size(percresp,2)),percresp);
% scatter(repmat(2,size(percnot,1),size(percnot,2)),percnot);
% 
% figure;
% hold all
% vecnm={'percresp','percnot'};
% lowbound=[];
% upbound=[];
% for inm=1:size(vecnm,2)
%     s=['[xb,pb]=empcdf(' vecnm{inm} ');'];
%     eval(s);
%     tmp=find(pb<0.95);
%     if isempty(tmp)
%         upbound(inm)=xb(max(find(pb==max(pb))));
%     else
%         upbound(inm)=xb(max(tmp));
%     end
%     tmp=find(pb<0.05);
%     if isempty(tmp)
%         lowbound(inm)=xb(max(find(pb==min(pb))));
%     else
%         lowbound(inm)=xb(max(tmp));
%     end
%     stairs(xb,pb,'LineWidth',3)
% end
% legend(vecnm);
% 
% figure;
% plotvec=[nanmean(percresp),nanmean(percnot)];
% errorbar([1,2],plotvec,plotvec-lowbound,upbound-plotvec,'ok');
% if kstest2(percresp,percnot)
%     line([1,2],[1.2,1.2],'color','r');
% end
% ylabel('percent low variability')
% set(gca,'YLim',[0,1.3],'XLim',[0,3])
% set(gca,'XTick',[1,2])
% set(gca,'XTickLabel',{'response','noresp'})
%%
db=unique(thisdb);
db=db(2:end);
for idb=1:size(db,2)
    theseinds=find(thisdb==db(idb));
%     figure;
   
%     hold all
    vecnm={'medbound(theseinds)','meanresp(theseinds)','meannot(theseinds)'};
    upbound=[];
    lowbound=[];
    for inm=1:size(vecnm,2)
        s=['[xb,pb]=empcdf(' vecnm{inm} ');'];
        eval(s);
        tmp=find(pb<0.85);
        if isempty(tmp)
            upbound(inm)=xb(max(find(pb==max(pb))));
        else
            upbound(inm)=xb(max(tmp));
        end
        tmp=find(pb<0.15);
        if isempty(tmp)
            lowbound(inm)=xb(max(find(pb==min(pb))));
        else
            lowbound(inm)=xb(max(tmp));
        end
%         stairs(xb,pb,'LineWidth',3)
    end
%     legend(vecnm);
%      title([num2str(db(idb)) ' dB SPL']);
    
    figure;
    
    hold on
    scatter(repmat(0,size(medbound(theseinds),1),size(medbound(theseinds),2)),medbound(theseinds));
    scatter(repmat(1,size(meanresp(theseinds),1),size(meanresp(theseinds),2)),meanresp(theseinds));
    scatter(repmat(2,size(meannot(theseinds),1),size(meannot(theseinds),2)),meannot(theseinds));
    title([num2str(db(idb)) ' dB SPL']);
    
    figure;
    
    plotvec=[nanmean(medbound(theseinds)),nanmean(meanresp(theseinds)),nanmean(meannot(theseinds))];
errorbar([0,1,2],plotvec,plotvec-lowbound,upbound-plotvec,'ok');
     if kstest2(medbound(theseinds),meanresp(theseinds))
        line([0,1],[1,1],'color','r');
    end
    if kstest2(medbound(theseinds),meannot(theseinds))
        line([0,2],[1.2,1.2],'color','r');
    end
    if kstest2(meannot(theseinds),meanresp(theseinds))
        line([1,2],[1.4,1.4],'color','r');
    end
    set(gca,'YLim',[0,1.5])
 set(gca,'XTick',[0,1,2])
set(gca,'XTickLabel',{'spont','response','noresp'})
title([num2str(db(idb)) ' dB SPL']);

% 
%     figure;
%     hold on
%     scatter(repmat(1,size(percresp(theseinds),1),size(percresp(theseinds),2)),percresp(theseinds));
%     scatter(repmat(2,size(percnot(theseinds),1),size(percnot(theseinds),2)),percnot(theseinds));
%       title([num2str(db(idb)) ' dB SPL']);
%     
%     
% figure;
% hold all
% vecnm={'percresp(theseinds)','percnot(theseinds)'};
% lowbound=[];
% upbound=[];
% for inm=1:size(vecnm,2)
%     s=['[xb,pb]=empcdf(' vecnm{inm} ');'];
%     eval(s);
%     tmp=find(pb<0.95);
%     if isempty(tmp)
%         upbound(inm)=xb(max(find(pb==max(pb))));
%     else
%         upbound(inm)=xb(max(tmp));
%     end
%     tmp=find(pb<0.15);
%     if isempty(tmp)
%         lowbound(inm)=xb(max(find(pb==min(pb))));
%     else
%         lowbound(inm)=xb(max(tmp));
%     end
%     stairs(xb,pb,'LineWidth',3)
% end
% legend(vecnm);
%  title([num2str(db(idb)) ' dB SPL']);
% 
% figure;
% plotvec=[nanmean(percresp(theseinds)),nanmean(percnot(theseinds))];
% errorbar([1,2],plotvec,plotvec-lowbound,upbound-plotvec,'ok');
% if kstest2(percresp(theseinds),percnot(theseinds))
%     line([1,2],[1.2,1.2],'color','r');
% end
% ylabel('percent low variability')
% set(gca,'YLim',[0,1.3],'XLim',[0,3])
% set(gca,'XTick',[1,2])
% set(gca,'XTickLabel',{'response','noresp'})
% title([num2str(db(idb)) ' dB SPL']);
end

%%
%before this will work... I need to put in some lines for normalizing each
%record to (the mean of baseline?) or to (the Vrest of baseline?... 
    %probably the later to show it is zero)
%plot baseline period in subplot to the left for each as well (from stimonset = 0) to compare
varfig = figure;
hold on
sigfig = figure;
hold on
tmp=[];
for ivar=1:size(allvardata,2)
    sigoff=allsigwin{ivar}(2);
    recover=round(0.7/expt.wc.dt);
    tmpvar(ivar,:)=allvardata{ivar}(1,sigoff:sigoff+recover);
    figure(varfig)
    line([1:size(tmpvar(ivar,:),2)]*expt.wc.dt,tmpvar(ivar,:)) 
     tmpsig(ivar,:)=allsigdata{ivar}(1,sigoff:sigoff+recover);
    figure(sigfig)
       line([1:size(tmpsig(ivar,:),2)]*expt.wc.dt,tmpsig(ivar,:)) 
end
figure(sigfig)
 plot(nanmean(tmpsig)')

%%
%to ask whether the neuron is in up/down states during baseline
...and if during response also or not
%get distribution of Vm for each cell across all stimuli for that cell
...test for binary (or poisson) distribution (whatever they use for spike widths)
    
for iexpt=1:size(repexpts,2)
    thisstim=allsignames{iexpt};
    thisexpt=repexpts{iexpt};
    load([r.Dir.Expt thisexpt])
    
    vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
    ...this is the potential all ic data recorded at in this set
        
sigexpt=filtesweeps(vmexpt,0,'wavnames',thisstim)
sigon=allsigwin{iexpt}(1);
sigoff=allsigwin{iexpt}(2);
tmpresp=sigexpt.wc.data(:,sigon:sigoff);
[xr,pr]=empcdf(tmpresp(:));

basetimes=allbasewin{iexpt};
tmpbase=sigexpt.wc.data(:,basetimes(1):basetimes(2));
[xb,pb]=empcdf(tmpbase(:));

figure;
hold on
stairs(xr,pr,'LineWidth',3,'color','r')
stairs(xb,pb,'LineWidth',3,'color','b')

edges=[-0.080:0.002:-0.030];
nbase=histc(tmpbase(:),edges)./size(tmpbase(:),1);
nresp=histc(tmpresp(:),edges)./size(tmpresp(:),1);
figure;
line(edges,nbase','color','b')
hold on
line(edges, nresp,'color','r')

end
%%

%% variables to clear for next section
trials = 8; %this way i can get the 80db stim from that one cell...
vmresp=[];
base_vmresp=[];
base_spkresp=[];
spkresp=[];
respDB=[];
respdata=[];
cb=[];
cr=[];
respslope=[];
baseslope=[];
%%
stimind=1;
respind=1;
hfig_vm = figure;
hold all
hfig_spk = figure;
hold all

use_dbinds = 0;
PlotResponseGrid = 0;

for iexpt=1:size(repexpts,2)
    iexpt
    thisexpt=repexpts{iexpt};
    load([r.Dir.Expt thisexpt])
    
    vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
    ...this is the potential all ic data recorded at in this set
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
        % get periods of Vm response across db
        up_win = [];
        low_win = [];
        confint = [];
        upinds = [];
        lowinds = [];
        up_inds = [];
        low_inds = [];
        spiketimes = [];
        spikesmat = [];
        nreps = [];
        b=[];
        sdata=[];
        conf_p = 95;
        windowsize = 500;
        binsize = 10;
        p = 0.85;
       
        %get the response and inhibited data off confint
        for idb=1:size(thisdb,2)
            sigexpt=filtesweeps(vmexpt,0,'wavnames',thiscond(idb).wavnames);
            sdata(idb,:)=mean(medfilt1(sigexpt.wc.data,200,[],2))*1000;
            ydatabound(idb,:) = [min(min(sdata(idb,basetimes(1):end))), max(max(sdata(idb,sigon:sigoff)))];
            b(idb)=mean(sdata(idb,basetimes(1):basetimes(2)));
            confint(idb,:) = getCDFconf (sdata(idb,basetimes(1):sigon),conf_p);
            
            [upinds{idb}, lowinds{idb}] = WindowResponse(sdata(idb,sigon:sigoff), confint(idb,:), windowsize, binsize, p);
            
            fs=1/expt.wc.dt;
            highpassdata=HighpassGeneral(sigexpt.wc.data,[],metadata.highcutoff,fs);
            if negative==1
                highpassdata=-highpassdata;
            end
            nreps(idb)=size(highpassdata,1);
            spikesmat{idb} = getspikesmat(highpassdata,spikesthresh,expt);
            for itrial=1:  nreps(idb)
                spiketimes{idb,itrial}=find(spikesmat{idb}(itrial,:));
            end
        end
        meanbase=mean(b);
        
        if use_dbinds ==1
            %can get union of response windows
            up_inds = union(up_inds, upinds);
            low_inds = union(low_inds, lowinds);
            %get windows from continuous indices
            if isempty(up_inds)
                up_win = [];
            elseif ~isempty(up_inds)
                up_win = getWindowEdges (up_inds, 1, 1);
                up_win=up_win+sigon;
            end
            
            if isempty(low_inds)
                low_win = [];
            elseif ~isempty(low_inds)
                low_win = getWindowEdges (low_inds, 1, 1);
                low_win=low_win+sigon;
            end
        end
        
        if use_dbinds == 0
            %can get windows based on average
            confint = getCDFconf (mean(sdata(:,basetimes(1):sigon),1),conf_p);
            [up_inds, low_inds] = ...
                WindowResponse(mean(sdata(:,sigon:sigoff),1), confint, ...
                windowsize, binsize, p);
            up_win = getWindowEdges (up_inds, 1, 1)+sigon;
            low_win = getWindowEdges (low_inds, 1, 1)+sigon;
        end
        
        if PlotResponseGrid ==1
            %plot the data and the windows
            respfig = figure;
            hold on
            scaleticks = 1;
            [grad,im]=colorGradient([0.08,0,0.5],[0,0.78,0.78],3);
            ymax = max(max(ydatabound));
            y_raster = ymax;
            ymin = min(min(ydatabound));
            xtime=[1:size(vmexpt.wc.data,2)]*expt.wc.dt;
            for idb=1:size(thisdb,2)
                line(xtime(1,basetimes(1):end),sdata(idb,basetimes(1):end),'color',grad(idb,:),'LineWidth',3);
                %             plot([xtime(1),xtime(end)],[confint(idb,2),confint(2)],'--','color','k')
                %add a line between each raster plot
                line([xtime(1),xtime(end)],[y_raster,y_raster],'color','k');
                text(xtime(1),y_raster,[num2str(thisdb(idb)) 'dB SPL'],...
                    'HorizontalAlignment','center',	'BackgroundColor', grad(idb,:),...
                    'color',[1,1,1]);
                
                for itrial=1:  nreps(idb)
                    for ispike=1:size(spiketimes{idb,itrial},2)
                        plot([spiketimes{idb,itrial}(ispike)*expt.wc.dt,...
                            spiketimes{idb,itrial}(ispike)*expt.wc.dt],...
                            [(itrial*scaleticks)+y_raster, ...
                            (scaleticks*(itrial+0.9))+y_raster], ...
                            'color','k','LineWidth',2)
                    end
                end
                y_raster = y_raster + nreps (idb) +1;
            end
            axis tight
            ylims = get(gca,'YLim');
            set(gca,'YLim',[ymin,ylims(2)],'YTick',...
                [(floor(ymin)-mod(floor(ymin),5)):5:(ceil(ymax)+mod(ceil(ymax),5))],...
                'XTick',[0:1:floor(xtime(end))],'TickDir','out')
            text(xtime(1),round(ymin),[num2str(round(ymin)) 'mV'],...
                'HorizontalAlignment','center',	'BackgroundColor', 'k',...
                'color',[1,1,1]);
            box off
            SigTimeBox(gca, sigon*expt.wc.dt,sigoff*expt.wc.dt, get(gca,'YLim'),[0.5 0.5 0.5]);
            SigTimeBox(gca, xtime(1),xtime(end), [mean(confint(:,1)),  mean(confint(:,2))],'k');
            for iresp=1:size(up_win,2)
                SigTimeBox(gca, up_win(1,iresp)*expt.wc.dt, ...
                    up_win(2,iresp)*expt.wc.dt, get(gca,'YLim'),'r');
            end
            for inhib=1:size(low_win,2)
                SigTimeBox(gca, low_win(1,inhib)*expt.wc.dt, ...
                    low_win(2,inhib)*expt.wc.dt, get(gca,'YLim'),'b');
            end
            title(['CELL: ' expt.name '   SIGNAL: ' thiscond(istimcond).wavnames(1:end-3)], 'Interpreter','none')
            set(respfig,'Position',[   252   531   957   404]);
            
            saveas(respfig,[r.Dir.Expt 'Analysis/Meta_10Trial_Responses_DB/' expt.name '_'...
                num2str(istimcond) '_Responses_Windows.fig']);
            saveas(respfig,[r.Dir.Expt 'Analysis/Meta_10Trial_Responses_DB/' expt.name '_'...
                num2str(istimcond) '_Responses_Windows.png']);
        end
        % need to get spiking responses that intersect with Vm responses...
        % not just look at spiking during Vm responses like i do here...
       
        normdata = sdata - meanbase;
        for iresp=1:size(up_win,2)
            trialspk = [];
            basespk = [];
            tmpvmresp = [];
            tmpbasevm = [];
            for idb=1:size(thisdb,2)
                tmpvmresp(idb,:)=normdata(idb,up_win(1,iresp):up_win(2,iresp));
                tmpbasevm(idb,:)=normdata(idb,basetimes(1):basetimes(2));
                
                for itrial=1:size(spiketimes,2)
                    highind = find(spiketimes{idb,itrial} > up_win(1,iresp));
                    lowind = find(spiketimes{idb,itrial} < up_win(2,iresp));
                    trialspk(idb,itrial) = size(intersect(highind,lowind)',1);
                    highind = find(spiketimes{idb,itrial} > basetimes(1));
                    lowind = find(spiketimes{idb,itrial} < basetimes(2));
                    basespk (idb,itrial) =  size(intersect(highind,lowind)',1);
                end
            end
            
            trialspk = trialspk ./ (diff(up_win(:,iresp)))*expt.wc.dt;
            trialspk = trialspk ./ max(max([trialspk;basespk]));
            nanind = isnan(trialspk);
            trialspk(nanind) = 0;
            basespk = basespk ./ diff(basetimes*expt.wc.dt);
            basespk = basespk ./ max(max([trialspk;basespk]));
            nanind = isnan(basespk);
            basespk(nanind) = 0;
            base_spkresp{respind} = mean(basespk,2);
            spkresp{respind}=nanmean(trialspk,2);
            cs{respind}=fit(thisdb',spkresp{respind},'poly1');
            spkslope(respind)=cs{respind}.p1;
            
            vmresp{respind}=mean(tmpvmresp./max(max(tmpvmresp)),2);
            respDB{respind}=thisdb;
            cr{respind}=fit(thisdb',vmresp{respind},'poly1');
            respslope(respind)=cr{respind}.p1;
            
            base_vmresp{respind} = mean(tmpbasevm./max(max(tmpbasevm)),2);
            cb{respind}=fit(thisdb',base_vmresp{respind},'poly1');
            baseslope(respind)=cb{respind}.p1;
         
            figure(hfig_vm)
            scatter(respDB{respind},vmresp{respind});
            figure(hfig_spk)
            scatter(respDB{respind},spkresp{respind});
            respind=respind+1;
        end
        stimind=stimind+1;
    end
    
end
%%

mined=min([baseslope,respslope]);
maxed=max([baseslope,respslope]);
edges=[mined:(maxed-mined)/20:maxed];
rbin=histc(respslope,edges)/size(respslope,2);
bbin=histc(baseslope,edges)/size(baseslope,2);
figure;
line(edges, rbin,'color','r')
line(edges, bbin,'color','b')

% how many responses increase with db?
% first find confidence bounds on baseline slopes
confint=getCDFconf(baseslope,85);
%find significantly positive response slopes
upperslope=find(respslope>confint(2));
lowerslope=find(respslope<confint(1));
[xr,pr]=empcdf(respslope);
[xb,pb]=empcdf(baseslope);
figure
hold on
stairs(xr,pr,'LineWidth',3,'color','r')
stairs(xb,pb,'LineWidth',3,'color','b')
title([num2str((size(upperslope,2)/size(respslope,2))*100) ' % inc slope ; ' ...
    num2str((size(lowerslope,2)/size(respslope,2))*100) ' % dec slope'] )


%%
figure;
hold all
respind=1;
tmp1=[];
tmp2=[];
tmp3=[];
bmp1=[];
bmp2=[];
bmp3=[];

dblist=[40,60,80];
for imdb=1:size(dblist,2)
    metadb=dblist(imdb);
    
    for irespind=1:size(vmresp,2)
        thisdb=respDB{irespind};
        for idb=1:size(thisdb,2)
            if (thisdb(idb)==metadb)
                dbresp{imdb}(irespind)=vmresp{irespind}(idb);
            end
        end
    end
end
allbase=[];
for ibase=1:size(base_vmresp,2)
  allbase=[allbase,base_vmresp{ibase}];  
end

figure;
hold on
cb=getCDFconf(allbase,85);
scatter([0,0],cb,'d','r');
scatter(0,mean(allbase),80,'k','fill');
for idb=1:size(dblist,2)
cb=getCDFconf(dbresp{idb},85);
scatter([dblist(idb),dblist(idb)],cb,'d','r');
scatter(dblist(idb),mean(dbresp{idb}),80,'k','fill');
end
set(gca,'YLim',[-0.05,1],'XLim',[-10,90],'XTick',[0,40,60,80],...
    'XTickLabel',{'baseline','40dB','60dB','80dB'})
ylabel('average fraction max Vm')
title([num2str(size(vmresp,2)) ' responses ; ' ...
    num2str(size(normbase,2)) ' signals ; ' ...
    num2str(size(repexpts,1)) ' cells'])
