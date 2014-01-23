function [sigdata,vardata,basewin] = GetVarDistrib(expt,stimcond,r)

% function [sigdata,vardata,onoff,percsignif,d] = GetVarDistrib(expt,stimcond,r)
% varbase{isig,iclamp}=var(vmdata(:,expt.analysis.params.baselinewin(1):sigon));
% plotmax(1)=max(varbase{isig,iclamp});
% varsig{isig,iclamp}=var(vmdata(:,sigon:sigoff));
% plotmax(2)=max(varsig{isig,iclamp});
% maxedge= max(plotmax)-rem(max(plotmax),histbin)+histbin;
% edges=[0:histbin:maxedge];
% [nbase]=histc(varbase{isig,iclamp},edges,2);
% nbase=nbase./max(nbase);
% [nsig]=histc(varsig{isig,iclamp},edges,2);
% nsig=nsig./max(nsig);
% figure;
% hold on
% line(edges,nbase,'color','b')
% line(edges,nsig,'color','r')
% hvardistIC(isig,iclamp)=kstest2(nsig,nbase);
%%
% iclamp=1;
% hfillplot=figure;
% hdistplot=figure;
% hvarplot=figure;
% for isig=1:size(stimcond,2)
%     sigon=round(expt.analysis.params.waveonset_time/expt.wc.dt);
%     siglength=round(max(size(stimcond(isig).wavs))/44100/expt.wc.dt);
%     sigoff=sigon+siglength;
%     sigexpt=filtesweeps(vmexpt,0,'wavnames',stimcond(isig).wavnames);
%     vmdata=medfilt1(sigexpt.wc.data,200,[],2).*1000; %converted to millivolts
%     xtime=[1:size(vmdata,2)]*expt.wc.dt;
%     stder=std(vmdata,1);
%     upperstd=mean(vmdata,1) + stder;
%     lowerstd=mean(vmdata,1) - stder;
%
%
%     figure(hfillplot);
%     hs=subplot(size(stimcond,2),1,isig);
%     plot(xtime,mean(vmdata,1)); %,'color',mycolors{i}S
%     [fillhandle,msg]=jbfill(xtime,upperstd,lowerstd,[0.5 0.5 0.5],'k',1,0.5);
%     axis tight
%     SigTimeBox(hs, sigon*expt.wc.dt, sigoff*expt.wc.dt, get(gca,'YLim'));
%     if isig~=size(stimcond,2)
%         set(gca,'XTick',[])
%     end
%
%     histbin=20;
%     %get variance distribution of baseline and signal
%     varbase{isig,iclamp}=var(vmdata(:,expt.analysis.params.baselinewin(1):sigon));
%     plotmax(1)=max(varbase{isig,iclamp});
%     varsig{isig,iclamp}=var(vmdata(:,sigon:sigoff));
%     plotmax(2)=max(varsig{isig,iclamp});
%     maxedge= max(plotmax)-rem(max(plotmax),histbin)+histbin;
%     edges=[0:histbin:maxedge];
%     [nbase]=histc(varbase{isig,iclamp},edges,2);
%     nbase=nbase./size(varbase{isig,iclamp},2);
%     [nsig]=histc(varsig{isig,iclamp},edges,2);
%     nsig=nsig./size(varsig{isig,iclamp},2);
%
%     hvardistIC(isig,iclamp)=kstest2(nsig,nbase);
%     figure(hvarplot);
%     hs=subplot(size(stimcond,2),1,isig);
%     hold on
%     line(edges,nbase,'color','b')
%     line(edges,nsig,'color','r')
%     if isig~=size(stimcond,2)
%         set(gca,'XTick',[])
%     end
%
%     edges=[-80:2:-30];
%     tmp=vmdata(:,expt.analysis.params.baselinewin(1):sigon);
%     tmp=tmp(:);
%     [nbase]=histc(tmp,edges);
%     nbase=nbase./size(tmp,1);
%     tmp=vmdata(:,sigon:sigoff);
%     tmp=tmp(:);
%     [nsig]=histc(tmp,edges);
%     nsig=nsig./size(tmp,1);
%     VmdistIC(isig,iclamp)=kstest2(nsig,nbase);
%
%
%     figure(hdistplot);
%     hs=subplot(size(stimcond,2),1,isig);
%     hold on
%     line(edges,nbase,'color','b')
%     line(edges,nsig,'color','r')
%     if isig~=size(stimcond,2)
%         set(gca,'XTick',[])
%     end
%      set(gca,'XLim',[-80 -30])
%
% end

%%

basewin=expt.analysis.params.baselinewin;
s=PlotAndAsk(expt.wc.data(:,basewin(1):basewin(2)),'basebegin')
cellfun(@eval,s);
if ~isempty(basebegin)
    basewin(1)=basewin(1)+basebegin;
end
% vardata=[];
% sigdata=[];
% offsetvar=[];
% onsetvar=[];
% onoff=[];
% psignifresp=[];
% psignifbase=[];
% persignif=[];
% [sigon,sigoff]=GetSigTimes(expt,stimcond,1);
% hfig=figure;hold on
% set(hfig,'Position',[  20        -436        1029        1228]);
% hfillplot=figure;hold on
% set(hfillplot,'Position',[  20        -436        1029        1228]);
% fcumm=figure;hold on
% set(fcumm,'Position',[ 20        -436        1029        1228]);
for isig=1:size(stimcond,2)
    sigexpt=filtesweeps(expt,0,'wavnames',stimcond(isig).wavnames);
    filtdata=medfilt1(sigexpt.wc.data,200,[],2);
%     
%     stder=std(filtdata,1);
%     upperstd=mean(filtdata,1) + stder;
%     lowerstd=mean(filtdata,1) - stder;
%     xtime=[1:size(filtdata,2)]*expt.wc.dt;
%     figure(hfillplot);
%     hs=subplot(size(stimcond,2),1,isig);
%     plot(xtime,mean(filtdata,1)); %,'color',mycolors{i}S
%     [fillhandle,msg]=jbfill(xtime,upperstd,lowerstd,[0.5 0.5 0.5],'k',1,0.5);
%     axis tight
%     SigTimeBox(hs, sigon*expt.wc.dt, sigoff*expt.wc.dt, get(gca,'YLim'));
%     if isig~=size(stimcond,2)
%         set(gca,'XTick',[])
%     end
%        title(stimcond(isig).wavnames,'Interpreter','none')
%        
    sigdata(isig,:)=mean(filtdata);
    vardata(isig,:)=var(filtdata);
    varallowt=round(0.1/expt.wc.dt); %use this filter to allow variance in real time to go above signif a little
%     vardataf=medfilt1(vardata(isig,:),varallowt,[],2);
vardataf=vardata(isig,:); % dont filter so much anymore... just start looking for offset after ~
%     basem=mean(vardataf(1,basewin(1):sigon));
%     basest=std(vardataf(1,basewin(1):sigon));
%     %     confbound=basem-basest; %~85% confbound
%     %*****cannot use standard dev as cutoff because not normal distrib!!!
%     [xb,pb]=empcdf(vardataf(1,basewin(1):sigon));
%     tmp=find(pb>0.15);
%     confbound=xb(min(tmp));
%     
%     figure(fcumm)
%     subplot(size(stimcond,2),1,isig);hold on
%     stairs(xb,pb,'color','b','LineWidth',3)
%     [xr,pr]=empcdf(vardataf(1,sigon:sigoff));
%     stairs(xr,pr,'color','r','LineWidth',3);
%     ylims=get(gca,'YLim');
%     line([confbound,confbound],ylims,'color','k','LineWidth',3)
%     title(stimcond(isig).wavnames,'Interpreter','none')
%    
%     d(isig)=kstest2(vardataf(1,basewin(1):sigon),vardataf(1,sigon:sigoff));
%     if d(isig)==0
%         onsetvar(isig)=nan;
%         offsetvar(isig)=nan;
%         psignifbase(isig)=nan;
%         psignifresp(isig)=nan;
%         
%         figure(hfig)
%         subplot(size(stimcond,2),1,isig);
%         hold on
%         line([basewin(1):size(vardataf,2)]*expt.wc.dt,vardataf(1,basewin(1):end)','color','k')
%         line([basewin(1),size(vardataf,2)]*expt.wc.dt,[confbound, confbound],'color',[0.5 0.5 0.5]);
%         axis tight
%         ylims=get(gca,'YLim');
%         SigTimeBox(gca,sigon*expt.wc.dt,sigoff*expt.wc.dt,ylims);
%         set(gca,'YTick',[])
%         if isig~=size(stimcond,2)
%             set(gca,'XTick',[])
%         end
%            title(stimcond(isig).wavnames,'Interpreter','none')
%     end
%     if d(isig)==1
%         figure(hfig)
%         subplot(size(stimcond,2),1,isig);
%         hold on
%         line([basewin(1):size(vardataf,2)]*expt.wc.dt,vardataf(1,basewin(1):end)','color','r')
%         line([basewin(1),size(vardataf,2)]*expt.wc.dt,[confbound, confbound],'color',[0.5 0.5 0.5]);
%         axis tight
%         ylims=get(gca,'YLim');
%         SigTimeBox(gca,sigon*expt.wc.dt,sigoff*expt.wc.dt,ylims);
%         set(gca,'YTick',[])
%         if isig~=size(stimcond,2)
%             set(gca,'XTick',[])
%         end
%         title(stimcond(isig).wavnames,'Interpreter','none')
%         
%         t=find(vardataf(1,basewin(1):sigon)<(confbound));
%         psignifbase(isig)=size(t,2)/size(vardataf(1,basewin(1):sigon),2);
%         t=find(vardataf(1,sigon:sigoff)<(confbound));
%         psignifresp(isig)=size(t,2)/size(vardataf(1,sigon:sigoff),2);
%         %downloaded crossing.m to find level crossings
%         ...http://www.mathworks.com/matlabcentral/fileexchange/2432-crossing
%             ind = crossing(vardataf,[],(confbound));
%         %gets the indices that the variance is lower than 1std (~85%) of the
%         %data; but only if the variance is at baseline level at signal onset
%         ontime=(ind(min(find(ind>sigoff)))-sigoff)*expt.wc.dt;
%         if ~isempty(ontime) && ontime~=0
%             onsetvar(isig)=ontime;
%         end
%         if isempty(ontime) || ontime==0
%             onsetvar(isig)=nan;
%         end
%         
%         %for now only get indices for when variance returns if the variance is
%         %low at the signal offset
%         %don't start looking for offset recovery until 100ms after stim
%         %offset
%         offtime=(ind(min(find(ind>sigoff+round(0.1/expt.wc.dt))))-sigoff)*expt.wc.dt;
%         if ~isempty(offtime) && offtime~=0
%             offsetvar(isig)=offtime;
%         end
%         %if the variance never recovers... record the max time recorded as
%         %the offset time
%         if isempty(offtime)
%            offsetvar(isig)=size(vardataf,2)*expt.wc.dt; 
%         end
%         if offtime==0
%             offsetvar(isig)=nan;
%         end
%         
%     end
end

%only use stimulus:response pairs in which the baseline variance is
... less than 15% + the standard deviation of signig percents during baseline for that cell
    % signifbnd=0.15+std(psignifbase);
% signifbnd=0.2;
% useinds=find(psignifbase<signifbnd);
%calculating differences and onset/offset, etc when d=1 for different distributions 
% useinds=find(d==1);
% 
% saveas(hfig,[r.Dir.Expt expt.name '\varPlot.fig'])
% saveas(hfig,[r.Dir.Expt expt.name '\varPlot.png'])
% 
%  saveas(fcumm,[r.Dir.Expt expt.name '\varCumDist.fig'])
%     saveas(fcumm,[r.Dir.Expt expt.name '\varCumDist.png'])
%     close(fcumm)
% 
% figure(hfillplot)
% set(gcf, 'Renderer', 'ZBuffer')
% saveas(hfillplot,[r.Dir.Expt expt.name '\fillPlot.fig'])
% saveas(hfillplot,[r.Dir.Expt expt.name '\fillPlot.png'])
% 
% close(hfillplot)
% close(hfig)
% 
% 
% hfig=figure;
% hold all
% mn=[nanmean(onsetvar(useinds)), nanmean(offsetvar(useinds))]*1000;
% st=[nanstd(onsetvar(useinds)),nanstd(offsetvar(useinds))]*1000;
% h(1)=bar(1,mn(1),'b');
% errorbar(1,mn(1),st(1),'color','k')
% h(2)=bar(2,mn(2),'r');
% errorbar(2,mn(2),st(2),'color','k')
% set(gca,'YLim',[0 max(mn)+max(st)],'XLim',[0 3],'XTick',[1 2],'XTickLabel',{'var onset' ,'var offset'})
% ylabel('milliseconds to decrease var after onset and inc after offset')
% title([expt.name ' ' num2str(size(useinds,2)) ' stimuli'],'Interpreter','none')
% legend(h,num2str(mn'))
% saveas(hfig,[r.Dir.Expt expt.name '\VarOnsetOffset.fig'])
% saveas(hfig,[r.Dir.Expt expt.name '\VarOnsetOffset.png'])
% close(hfig)
% onoff=[onsetvar', offsetvar'];
% 
% hfig=figure;
% hold all
% st=[std(psignifbase(useinds)),std(psignifresp(useinds))];
% mn=[mean(psignifbase(useinds)),mean(psignifresp(useinds))];
% h(1)=bar(1,mn(1),'b');
% errorbar(1,mn(1),st(1),'color','k')
% h(2)=bar(2,mn(2),'r');
% errorbar(2,mn(2),st(1),'color','k')
% set(gca,'YLim',[0 1+max(st)],'XLim',[0 3],'XTick',[1 2],'XTickLabel',{'baseline' ,'stimulus'})
% ylabel('fraction variance below 1SD significance level')
% title([expt.name ' ' num2str(size(useinds,2)) ' stimuli'],'Interpreter','none')
% legend(h,num2str(mn'))
% saveas(hfig,[r.Dir.Expt expt.name '\FracSignifVar.fig'])
% saveas(hfig,[r.Dir.Expt expt.name '\FracSignifVar.png'])
% close(hfig)
% percsignif=[psignifbase', psignifresp'];

