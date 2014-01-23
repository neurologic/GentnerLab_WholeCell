function subplotMotifByVm(expt,Vm,ncol,lowpass,params,havestimcond,nrepsstim,folder,plotrasterinfo,plotrestricted,Yunits,r,saveme,meanonly)
% plotraster=0;
% plotrestricted=1;
% aftermotif=0.75;
% sortstims.fieldnames='stimtype';
% sortstims.sortvalues={''};

waveonset_time = params.waveonset_time;
baselinewin=params.baselinewin;

isplotraster=plotrasterinfo.isplotraster;
if isplotraster==1
    spikesVm=plotrasterinfo.spikesVm;
    %spikesthresh=plotrasterinfo.spikesthresh;
    cutoff=plotrasterinfo.cutoff;
    %scaleticks=plotrasterinfo.scaleticks;
    %negative=plotrasterinfo.negative;
    spikesexpt=filtesweeps(expt,0,'Vm',spikesVm);
end

isplotrestricted=plotrestricted.isplotrestricted;
if isplotrestricted==1
    AfterMotDur=plotrestricted.AfterMotDur;
end


clear mots %cond
clear Vmexpt
clear h

cond.Vm.cond.fn = @(sweeps) sweeps.Vm==Vm;
cond.Vm.cond.Vm = Vm;
cond.Vm.cond.color = 'k';
cond.Vm.cond.name = 'PSC';
cond.Vm.cond.range =  baselinewin(1);%1; %

Vmexpt = filtesweeps(expt,0,'Vm',cond.Vm.cond.Vm);
% Vmexpt = addwaves(Vmexpt,waveonset_time,folder,[]);

if isplotraster==1
    highpassdata=HighpassGeneral(spikesexpt.wc.data,baselinewin,[],cutoff,1/expt.wc.dt);
    figure;
    hold on;
    plot(highpassdata');
    spikesthresh=input('threshold');
    negative=input('invert spikes? 1 or 0');
    scaleticks=spikesthresh;
end
if isempty(havestimcond)
    sortstims.fieldnames='Vm';
    sortstims.sortvalues=Vm;
    thisstimcond = filtstimcond(Vmexpt,folder,waveonset_time,nrepsstim,sortstims);
end
thiswavlist=unique(Vmexpt.sweeps.wavnames);
thisstimcond=getsubstimcond(havestimcond,thiswavlist);

nrow = ceil(length(thisstimcond)/ncol);
xtime = single([1:size(Vmexpt.wc.data,2)]*Vmexpt.wc.dt);

h.hfig = figure('Visible','off');
hold on
for imotif = 1:length(thisstimcond)
    
    h.hax(imotif) = subplot(nrow,ncol,imotif);
    hold on
    
    this_expt = filtesweeps(Vmexpt,0,'wavnames',thisstimcond(imotif).wavnames);
    w = thisstimcond(imotif).wavs;
    
    trials = find(cond.Vm.cond.fn(this_expt.sweeps)); 
%    this_data = this_expt.wc.data(trials,:); dont know what these
%     steps were for...
    this_data = this_expt.wc.data;
    
    %this_data = BaselineGeneral(this_data,baselinewin,[]);
    if lowpass==1
        this_data = LowpassGeneral(this_data,params.baselinewin,[]);
    end
    
    xrange=cond.Vm.cond.range;
    
    if isplotrestricted==1
        sweepstop = round((waveonset_time+AfterMotDur)/expt.wc.dt)+round(length(w)/44100/expt.wc.dt);
        %hl =
        if meanonly==0
            h2 = line(xtime(xrange:sweepstop),this_data(:,xrange:sweepstop)','color','b');%cond.Vm.cond(iVm).color
             hl = line(xtime(xrange:sweepstop),mean(this_data((1:size(trials,1)),xrange:sweepstop)),'color','k','LineWidth',1);
        else
            hl = line(xtime(xrange:sweepstop),mean(this_data((1:size(trials,1)),xrange:sweepstop)),'color','k','LineWidth',1);
        end
    else
        if meanonly==0
            h2 = line(xtime(xrange:end),this_data(:,xrange:end)','color','b');%cond.Vm.cond(iVm).color
            hl = line(xtime(xrange:end),mean(this_data((1:size(trials,1)),xrange:end)),'color','k','LineWidth',1);
        else
            hl = line(xtime(xrange:end),mean(this_data((1:size(trials,1)),xrange:end)),'color','k','LineWidth',1);%[0.5 0.5 0.5]
        end
    end
    axis tight
    ymaxoffset=max(max(this_data(:,xrange:end),[],1));
    yminoffset=min(min(this_data(:,xrange:end),[],1));
    if isplotraster==1
        plotraster_generalExpt(h.hax(imotif),imotif,ymaxoffset,scaleticks,expt,spikesVm,negative,cutoff,spikesthresh,baselinewin,thisstimcond);
    end
    axis tight
    %         if imotif==1
    %             set(gca,'TickDir','out');
    %         else
    %set(gca,'YTick',[],'YTickLabel',[]);
    set(gca,'TickDir','out','XTick',[],'XTickLabel',[]); %,'YTick',[],'YTickLabel',[]
    %         end
    axis tight
    if isplotraster==1
        addstimtoplot('bottom',imotif,h.hax(imotif),thisstimcond,params,Vmexpt,1);
    else
        addstimtoplot('top',imotif,h.hax(imotif),thisstimcond,params,Vmexpt,1);
    end
%     title(thisstimcond(imotif).wavnames{:});
    title(thisstimcond(imotif).wavnames);

axis tight
    xlim=get(gca,'XLim');
  
    set(gca,'XLim',[xtime(xrange) xlim(2)]);
%     set(gca,'YLim',[-100 30]);
end
set(h.hfig ,'Visible','on')
set(h.hfig,'Position',[[116 460 1019 898]]);%[318 1149 1053 801]
%     setXLabel(h.hax(1),'sec');
%     setYLabel(h.hax(1),Yunits);
if saveme==1
    if lowpass==1
        saveas(h.hfig,fullfile(r.Dir.Expt,expt.name,[expt.name(4:end) '_' num2str(Vm) 'lowpass.fig']))
        saveas(h.hfig,fullfile(r.Dir.Expt,expt.name,[expt.name(4:end) '_' num2str(Vm) 'lowpass.jpg']))
        saveas(h.hfig,fullfile(r.Dir.Expt,expt.name,[expt.name(4:end) '_' num2str(Vm) 'lowpass.eps']))
        if meanonly==1
            saveas(h.hfig,fullfile(r.Dir.Expt,expt.name,[expt.name(4:end) '_' num2str(Vm) 'lowpassMeanonly.fig']))
        saveas(h.hfig,fullfile(r.Dir.Expt,expt.name,[expt.name(4:end) '_' num2str(Vm) 'lowpassMeanonly.jpg']))
        saveas(h.hfig,fullfile(r.Dir.Expt,expt.name,[expt.name(4:end) '_' num2str(Vm) 'lowpassMeanonly.eps']))
        end
    else if isplotraster==1
            saveas(h.hfig,fullfile(r.Dir.Expt,expt.name,[expt.name(4:end) '_' num2str(Vm) 'raster.fig']))
            saveas(h.hfig,fullfile(r.Dir.Expt,expt.name,[expt.name(4:end) '_' num2str(Vm) 'raster.jpg']))
            saveas(h.hfig,fullfile(r.Dir.Expt,expt.name,[expt.name(4:end) '_' num2str(Vm) 'raster.eps']))
        else
            saveas(h.hfig,fullfile(r.Dir.Expt,expt.name,[expt.name(4:end) '_' num2str(Vm) '.fig']))
            saveas(h.hfig,fullfile(r.Dir.Expt,expt.name,[expt.name(4:end) '_' num2str(Vm) '.jpg']))
            saveas(h.hfig,fullfile(r.Dir.Expt,expt.name,[expt.name(4:end) '_' num2str(Vm) '.eps']))
        end
    end
end