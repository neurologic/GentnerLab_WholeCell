r = rigdef('mac');
% load(fullfile([r.Dir.Expt 'KP_B130_131120_p1c2a.mat']));
load(fullfile([r.Dir.Expt 'KP_B694_131212_p1c1.mat']));

gsyn_B694_p1c1 = gsyn;
esyn_B694_p1c1 = esyn;
ci_B694_p1c1 = ci;

load('/Users/kperks/GitHub/Data_Mat/Analysis/G_E_Anal/G_E_Anal.mat')
%% calculate gsyn and esyn from VC holding clamps without correcting for baseline
VC=[-40,-60,-80];
clearvars gsyn esyn ci
[grad,im]=colorGradient([1,0,0],[0,0,1],size(VC,2));
[dbstimcond,dblevels]=getDBstimcond(expt);

stimcond=dbstimcond{1};
% stimcond = expt.stimcond;
for isig=1:size(expt.stimcond,2)
    [sigon,sigoff]=GetSigTimes(expt,stimcond,isig);
    sigexpt=filtesweeps(expt,0,'wavnames',stimcond(isig).wavnames);
    
    thisresp=[];
    holding = [];
    vmdata = [];
    for iclamp=1:max(size(VC))
        vmexpt=filtesweeps(sigexpt,0,'Vm',VC(iclamp));
        vmdata{iclamp}=medfilt1(vmexpt.wc.data,200,[],2);
        holding{iclamp} = repmat(VC(iclamp),size(vmdata{iclamp},1),size(vmdata{iclamp},2));
    end

    for isamp = 1:size(vmdata{1},2)
        theseclamps = [];
        thesedata = [];
        for iclamp=1:max(size(VC))
            theseclamps = [theseclamps; holding{iclamp}(:,isamp)];
           thesedata = [thesedata; vmdata{iclamp}(:,isamp)];
        end
         cf=fit(theseclamps,thesedata,'poly1');
            ci{isig,isamp}=round(confint(cf));
            gsyn(isig,isamp)= cf.p1;
            esyn(isig,isamp)=-cf.p2/cf.p1;
    end
    
    
    
end
for s = 1:size(ci,1)
    for i = 1:size(ci,2)
        gsyn_ci_up(s,i) = ci{s,i}(1,1);
        gsyn_ci_dn(s,i) = ci{s,i}(2,1);
        esyn_ci_up(s,i) = ci{s,i}(1,2)./gsyn(s,i);
        esyn_ci_dn(s,i) = ci{s,i}(2,2)./gsyn(s,i);
    end
end
%%
save('/Users/kperks/GitHub/Data_Mat/Analysis/G_E_Anal/G_E_Anal.mat')

%%
%after filtering expt for only one set of signals...
% get respwindows from [meanresp, stdresp, respsig, respwindows, DBind]=SpikingResponseAcrossDB(expt,0)
VC=[-40,-60,-80];

% can get baseline manually...
medbase=[];
for iclamp=1:max(size(VC))
    vmexpt=filtesweeps(expt,0,'Vm',VC(iclamp));
    tmpdata=vmexpt.wc.data(:,vmexpt.analysis.params.baselinewin(1):vmexpt.analysis.params.baselinewin(2));
   tmpdata=medfilt1(tmpdata,200,[],2);
    tmp=[];
    for itmp=1:size(tmpdata,1)
       tmp=[tmp,tmpdata(itmp,:)];
    end
%     tmpdata=medfilt1(tmp,200,[],2);
    [varargout]=PlotAndAsk(tmp,'baseline');
    cellfun(@eval,varargout);
    medbase(iclamp)=baseline;
end

gsyn=[];
 esyn=[];
hfig=figure;hold on
 
for irespwin=1:size(respwindows,1)
    thiswint=diff(respwindows(irespwin,:));
    thiswin=round(respwindows(irespwin,:)/expt.wc.dt);
    
    for isig=1:size(stimcond,2)
        subplot(size(stimcond,2),1,isig);
        sigon=round(expt.analysis.params.waveonset_time/expt.wc.dt);
        siglength=round(max(size(stimcond(isig).wavs))/44100/expt.wc.dt);
        sigoff=sigon+siglength;
        sigexpt=filtesweeps(expt,0,'wavnames',stimcond(isig).wavnames);
        
        thisresp=[];
        [grad,im]=colorGradient([1,0,0],[0,0,1],size(VC,2));
       
        for iclamp=1:max(size(VC))
            vmexpt=filtesweeps(sigexpt,0,'Vm',VC(iclamp));
            %         vmdata=vmexpt.wc.data;
            vmdata=medfilt1(vmexpt.wc.data,200,[],2)-medbase(iclamp);
            thisresp(iclamp)=round(sum(mean(vmdata(:,thiswin(1):thiswin(2))),2));
            line([thiswin(1):thiswin(2)]*expt.wc.dt,mean(vmdata(:,thiswin(1):thiswin(2)))','color',(grad(iclamp,:)),'LineWidth',3)
        end
        legend(num2str(VC'))
       
        
        cf=fit(VC',thisresp','poly1');
        ci=confint(cf);
        ci=round(ci);
%         if size(find(ci(:,1)>0),1)==2
            gsyn(irespwin,isig)= cf.p1;
            esyn(irespwin,isig)=-cf.p2/cf.p1;
%         end
%         if size(find(ci(:,1)>0),1)~=2
%             gsyn(iresp,isig)= nan;
%             esyn(iresp,isig)=nan;
%         end
    end
    
end
[grad,im]=colorGradient([1,0,0],[0,0,1],size(respwindows,1));
 figure
 hold on
for iresp=1:size(respwindows,1)
   
scatter(dblevels{1},gsyn(iresp,:)/max(gsyn(iresp,:)),'MarkerFaceColor',grad(iresp,:),'MarkerEdgeColor','k')
set(gca,'YLim',[-0.01 1.01],'XLim',[35 85])
ylabel('percent max')
xlabel('dB SPL')
end
legend(num2str([1:size(respwindows,1)]))

   figure
 hold on
for iresp=1:size(respwindows,1)
 
scatter(dblevels{1},esyn(iresp,:),'MarkerFaceColor',grad(iresp,:),'MarkerEdgeColor','k')
set(gca,'XLim',[35 85])
ylabel('Esyn')
xlabel('dB SPL')
end
legend(num2str([1:size(respwindows,1)]'))

  [grad,im]=colorGradient([0,0,0.6],[0,1,1],3);
hfig=figure;hold on
vmexpt=filtesweeps(expt,0,'Vm',0)
for irespwin=1:size(respwindows,1)
    thiswint=diff(respwindows(irespwin,:));
    thiswin=round(respwindows(irespwin,:)/expt.wc.dt);
    for isig=1:size(stimcond,2)
       
        sigon=round(vmexpt.analysis.params.waveonset_time/expt.wc.dt);
        siglength=round(max(size(stimcond(isig).wavs))/44100/expt.wc.dt);
        sigoff=sigon+siglength;
        sigexpt=filtesweeps(vmexpt,0,'wavnames',stimcond(isig).wavnames);
        data=medfilt1(sigexpt.wc.data,200,[],2);
        line([thiswin(1):thiswin(2)]*expt.wc.dt,mean(data(:,thiswin(1):thiswin(2)))','color',(grad(isig,:)),'LineWidth',3)
        
    end
end
%%

[dbstimcond,dblevels]=getDBstimcond(expt);

stimcond=dbstimcond{1};
avgbase=[];
avgresp=[];
refVm=[];
VC=[-40,-60,-80,-100];%,-40];


% can get baseline manually...
medbase=[];
for iclamp=1:max(size(VC))
    vmexpt=filtesweeps(expt,0,'Vm',VC(iclamp));
    tmpdata=vmexpt.wc.data(:,vmexpt.analysis.params.baselinewin(1):vmexpt.analysis.params.baselinewin(2));
   tmp=[];
    for itmp=1:size(tmpdata,1)
       tmp=[tmp,tmpdata(itmp,:)];
    end
%     tmpdata=medfilt1(tmp,200,[],2);
    [varargout]=PlotAndAsk(tmp,'baseline');
    cellfun(@eval,varargout);
    medbase(iclamp)=baseline;
end

% and/or can get baseline by estimating based on Rin


%need to make sure that the RS are the same under the conditions....


avgbase=[];
avgresp=[];
avgbaseresp=[];
respslope=[];
respVCint=[];
hvardistVC=[];


for isig=1:size(stimcond,2)
    sigon=round(expt.analysis.params.waveonset_time/expt.wc.dt);
    siglength=round(max(size(stimcond(isig).wavs))/44100/expt.wc.dt);
    sigoff=sigon+siglength;
    sigexpt=filtesweeps(expt,0,'wavnames',stimcond(isig).wavnames);
    
    for iclamp=1:max(size(VC))
        vmexpt=filtesweeps(sigexpt,0,'Vm',VC(iclamp));
%         vmdata=vmexpt.wc.data;
        vmdata=medfilt1(vmexpt.wc.data,200,[],2);
        %need to use actual avg base waveforms so that sample by sample i
        %can check the slope bseteen the points... to get a distribution of
        %slopes to test the response slopes against...
        
        avgbase{isig}(iclamp,:)=max(mean(vmdata(:,expt.analysis.params.baselinewin(1):sigon)));
        avgbaseresp{isig}(iclamp,:)=mean(vmdata(:,expt.analysis.params.baselinewin(1):sigon))-medbase(iclamp);%-avgbase{isig}(iclamp,:);
        avgresp{isig}(iclamp,:)=mean(vmdata(:,sigon:sigoff),1)-medbase(iclamp);%-avgbase{isig}(iclamp,:);
%         
%         %get variance distribution of baseline and signal
%         varbase{isig,iclamp}=std(vmdata(:,expt.analysis.params.baselinewin(1):sigon));
%         plotmax(1)=max(varbase{isig,iclamp});
%         varsig{isig,iclamp}=std(vmdata(:,sigon:sigoff));
%         plotmax(2)=max(varsig{isig,iclamp});
%         maxedge= max(plotmax)-rem(max(plotmax),10)+10;
%         edges=[0:10:maxedge];
%         [nbase]=histc(varbase{isig,iclamp},edges,2);
%         nbase=nbase./max(nbase);
%         [nsig]=histc(varsig{isig,iclamp},edges,2);
%         nsig=nsig./max(nsig);
%         %these distributions can be plotted if needed
%         %      figure;
%         %      hold on
%         %      line(edges,nbase,'color','b')
%         %      line(edges,nsig,'color','r')
%         hvardistVC(isig,iclamp)=kstest2(nsig,nbase);
        
    end
    %get the distribution of slopes between holding potentials for each
    %baseline point
    for ibasesamp=1:size(avgbaseresp{isig},2)
        cf=fit(VC',avgbaseresp{isig}(:,ibasesamp),'poly1');
        baseslope(ibasesamp)=cf.p1;
    end
    baseedges=[(-(rem(min(baseslope),0.5)-min(baseslope)+0.5)):0.5:(max(baseslope)-rem(max(baseslope),0.5)+0.5)];
    nbaseslope=histc(baseslope,baseedges);
    [xb,pb]=empcdf(baseslope);
    tmp=find(pb>0.98);
    confbound=xb(min(tmp));
    figure;
    hold on
    stairs(xb,pb,'color','b','LineWidth',3)
    
    signifslope=[];
    %test if each slope during response is in or our of this distribution
    for irespsamp=1:size(avgresp{isig},2)
        cf=fit(VC',avgresp{isig}(:,irespsamp),'poly1');
        respslope{isig}(1,irespsamp)=cf.p1;
        respVCint{isig}(1,irespsamp)=-cf.p2/cf.p1;
        %test of nulll....
        if cf.p1>confbound
            signifslope(irespsamp)=1;
        else
            signifslope(irespsamp)=NaN;
        end
    end
    respedges=[(-(rem(min(respslope{isig}),0.5)-min(respslope{isig})+0.5)):0.5:(max(respslope{isig})-rem(max(respslope{isig}),0.5)+0.5)];
    nrespslope=histc(respslope{isig},respedges);
    [xr,pr]=empcdf(respslope{isig});
    stairs(xr,pr,'color','r','LineWidth',3);
    
    %          these distributions can be plotted if needed
    figure;
    hold on
    line(baseedges,nbaseslope./max(nbaseslope),'color','b','LineWidth',3)
    line(respedges,nrespslope./max(nrespslope),'color','r','LineWidth',3)
    legend('baseline "G"','response "G"')
    
    signifinds=find(isnan(signifslope)==0);
    notsignifinds=find(isnan(signifslope)==1);
    
    xtime=[1:size(respVCint{isig},2)]*expt.wc.dt;
    respVCint{isig}(1,notsignifinds)=NaN;
    respslope{isig}(1,notsignifinds)=NaN;
    
    %will want to compare the variance during synaptic input to the
    %variance not during synaptic input... and compare both to
    %baseline?
    %will i want to do this in the current clamp trace? if so, then do
    %i need to know how long after the synaptic input that the neuron
    %"responds"? like what is the time constant of the membrane...
    
end
xtime=[1:size(respVCint{isig},2)]*expt.wc.dt;
[grad,im]=colorGradient([0,0,0.6],[0,1,1],3);
figure;
hold on
for isig=1:size(stimcond,2)
    line(xtime,respVCint{isig},'color',(grad(isig,:)),'LineWidth',3);
end
legend('40','60','80dB')
figure;
hold on
for isig=1:size(stimcond,2)
    line(xtime,respslope{isig},'color',(grad(isig,:)),'LineWidth',3);
end
legend('40','60','80dB')

[grad,im]=colorGradient([1,0,0],[0,0,1],size(avgresp{isig},1));
for isig=1:size(avgresp,2)
    figure;
    hold all
    for iclamp=1:size(avgresp{isig},1)
        line([1:size(avgresp{isig},2)]*expt.wc.dt,avgresp{isig}(iclamp,:),'color',(grad(iclamp,:)),'LineWidth',3)
    end
    legend(num2str(VC'))
    title(stimcond(isig).wavnames);
end

for isig=1:size( avgbaseresp,2)
    figure;
    hold all
    for iclamp=1:size( avgbaseresp{isig},1)
        line([1:size( avgbaseresp{isig},2)]*expt.wc.dt, avgbaseresp{isig}(iclamp,:),'color',(grad(iclamp,:)),'LineWidth',4)
    end
    legend(num2str(VC'))
    title(stimcond(isig).wavnames);
end




%%
IC=[0];
medbaseIC=[];
avgrespIC=[];
avgbasestdIC=[];
sigrespindsIC=[];
notindsIC=[];
signifvarIC=[];
hvardistIC=[];
% stimcond=expt.stimcond;
for iclamp=1:max(size(IC))
    
    
    vmexpt=filtesweeps(expt,0,'Vm',IC(iclamp));
    hfillplot=figure;
    
    hold on
    title([num2str(IC(iclamp)) ' ' stimcond(1).wavnames])
    
    %find spike trheshold etc
    %cutoff frequency for highpass is 100
    
    
    for isig=1:size(stimcond,2)
        sigon=round(expt.analysis.params.waveonset_time/expt.wc.dt);
        siglength=round(max(size(stimcond(isig).wavs))/44100/expt.wc.dt);
        sigoff=sigon+siglength;
        sigexpt=filtesweeps(vmexpt,0,'wavnames',stimcond(isig).wavnames);
        vmdata=medfilt1(sigexpt.wc.data,200,[],2).*1000;
        xtime=[1:size(vmdata,2)]*expt.wc.dt;
        stder=std(vmdata,1);
        upperstd=mean(vmdata,1) + stder;
        lowerstd=mean(vmdata,1) - stder;
        
        
        figure(hfillplot);
        hs=subplot(size(stimcond,2),1,isig);
        plot(xtime,mean(vmdata,1)); %,'color',mycolors{i}S
        [fillhandle,msg]=jbfill(xtime,upperstd,lowerstd,[0.5 0.5 0.5],'k',1,0.5);
        axis tight
        SigTimeBox(hs, sigon*expt.wc.dt, sigoff*expt.wc.dt, get(gca,'YLim'));
        
        histbin=5;
        %get variance distribution of baseline and signal
        varbase{isig,iclamp}=var(vmdata(:,expt.analysis.params.baselinewin(1):sigon));
        plotmax(1)=max(varbase{isig,iclamp});
        varsig{isig,iclamp}=var(vmdata(:,sigon:sigoff));
        plotmax(2)=max(varsig{isig,iclamp});
        maxedge= max(plotmax)-rem(max(plotmax),histbin)+histbin;
        edges=[0:histbin:maxedge];
        [nbase]=histc(varbase{isig,iclamp},edges,2);
        nbase=nbase./size(varbase{isig,iclamp},2);
        [nsig]=histc(varsig{isig,iclamp},edges,2);
        nsig=nsig./size(varsig{isig,iclamp},2);
                figure;
                hold on
                line(edges,nbase,'color','b')
                line(edges,nsig,'color','r')
        hvardistIC(isig,iclamp)=kstest2(nsig,nbase);
        
        
        %get cummulative plot to find %98 bound on baseline
        [xb,pb]=empcdf(varbase{isig,iclamp});
        tmp=find(pb<0.05);
        confbound=xb(min(tmp));
                figure;
                hold on
                stairs(xb,pb,'color','b','LineWidth',3)
%         use confbound to see where variance is less than during baseline
        for irespsamp=1:size(varsig{isig,iclamp},2)
            
            if varsig{isig,iclamp}(1,irespsamp)<confbound
                signifvarIC{isig,iclamp}(1,irespsamp)=1;
            else
                signifvarIC{isig,iclamp}(1,irespsamp)=NaN;
            end
        end
        
        
        [xr,pr]=empcdf(varsig{isig,iclamp});
                stairs(xr,pr,'color','r','LineWidth',3);
        
        respdatavec=reshape(vmdata(:,expt.analysis.params.baselinewin(1):sigon),1,...
            size(vmdata,1)*size(vmdata(:,expt.analysis.params.baselinewin(1):sigon),2));
        [xbasevals,pbasevals]=empcdf(respdatavec );
        basedatavec=reshape( vmdata(:,sigon:sigoff),1,...
            size(vmdata,1)*size(vmdata(:,sigon:sigoff),2));
        [xrespvals,prespvals]=empcdf(basedatavec);
        %          these distributions can be plotted if needed
        
        medbaseIC{isig}(iclamp,:)=median(mean(vmdata(:,expt.analysis.params.baselinewin(1):sigon)));
        avgbasestdIC{isig}(iclamp,:)=mean(std(vmdata(:,expt.analysis.params.baselinewin(1):sigon)));
        
        
        avgrespIC{isig}(iclamp,:)=mean(vmdata(:,sigon:sigoff),1);%-medbase{isig}(iclamp,:);
        resptmp=avgrespIC{isig}(iclamp,:);
        medtmp=medbaseIC{isig}(iclamp,:);
        stdtmp=avgbasestdIC{isig}(iclamp,:);
        sigrespindsIC{isig,iclamp}= find(resptmp>(medtmp+(2*stdtmp)));
        notindsIC{isig,iclamp}=find(resptmp<(medtmp+(2*stdtmp)));
        
    end
    [grad,im]=colorGradient([0,0,0.6],[0,1,1],3);
    figure;
    hold all
    for isig=1:size(stimcond,2)
        line([1:size(avgrespIC{isig}(iclamp,:),2)]*expt.wc.dt,avgrespIC{isig}(iclamp,:),'color',(grad(isig,:)),'LineWidth',4)
    end
    % legend(num2str(clampval))
    title(num2str(IC(iclamp)))
    legend('40DB','60','80');
    
    figure;
    hold all
    for isig=1:size(stimcond,2)
        sigrespsIC=avgrespIC{isig}(iclamp,:);
        sigrespsIC(notindsIC{isig,iclamp})=NaN;
        line([1:size(sigrespsIC,2)]*expt.wc.dt,sigrespsIC,'color',(grad(isig,:)),'LineWidth',4)
    end
    %
    %  figure;
    % hold all
    % for isig=1:size(stimcond,2)
    %     VCintsub=VCint(isig,:);
    %     VCintsub(notindsIC{isig}(iclamp,:))=NaN;
    %
    % line([1:size(VCintsub,2)]/expt.wc.dt,VCintsub,'color',(grad(isig,:)),'LineWidth',4)
    % end
    %
    %  figure;
    % hold all
    % for isig=1:size(stimcond,2)
    %     slopesub=slope(isig,:);
    %     slopesub(notindsIC{isig}(iclamp,:))=NaN;
    %
    % line([1:size(slopesub,2)]/expt.wc.dt,slopesub,'color',(grad(isig,:)),'LineWidth',4)
    % end
    %
end
%
%
% for isig=1:size(stimcond,2)
%
% figure;
% hold all
% for iclamp=1:size(avgrespIC{isig},1)
%     sigrespsIC=avgrespIC{isig}(iclamp,:);
% line([1:size(sigrespsIC,2)]*expt.wc.dt,sigrespsIC,'color',(grad(iclamp,:)),'LineWidth',4)
% end
% end
