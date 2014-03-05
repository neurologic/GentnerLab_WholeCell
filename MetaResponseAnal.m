r=rigdef('mac')
%%


save([r.Dir.Expt 'Analysis/Meta_Responses_DB/exptsForICanal.mat'])
% ,...
%     'exptnames','readyexpt','sigdata','vardata','basewin','allstims','repexpts',...
%     'addstimcond','vmresp','base_vmresp','base_spkresp','spkresp','respDB',...
%     'respdata','cb','cr','respslope','baseslope','nresid_stim','nresid_base',...
%     's_base','s_stim','base_spk_vm','stim_spk_vm','stim_spk_thresh','stim_spk_thresh')

%%
load([r.Dir.Expt 'Analysis/Meta_Responses_DB/exptsForICanal.mat'])

%%
repexpts = {
    'KP_B130_131120_p1c1a.mat'
    'KP_B130_131120_p1c2a.mat'
    'KP_B130_131120_p1c2b.mat'
    'KP_B130_131120_p1c2c.mat'
    'KP_B130_131120_p1c3.mat'
    'KP_B136_131205_p1c2.mat'
    'KP_B136_131205_p2c2.mat'
    'KP_B136_131205_p3c1.mat'
    'KP_B580_120824_p1c2_b.mat'
    'KP_B680_131219_p1c1.mat'
    'KP_B680_131219_p1c2.mat'
    'KP_B682_130219_p1c1.mat'
    'KP_B689_131407_p2c1b.mat'
    'KP_B689_131407_p2c1c.mat'
    'KP_B689_131407_p2c1d.mat'
    'KP_B689_131407_p3c1.mat'
    'KP_B694_131212_p1c1.mat'
    'KP_B694_131212_p1c1b.mat'
    'KP_B694_131212_p2c1.mat'
    'KP_B694_131212_p2c2.mat'
    'KP_B694_131212_p3c1.mat'
    'KP_B694_131212_p4c1.mat'
    'KP_B790_140127_p1c1.mat'
    'KP_B790_140127_p1c2.mat'
    'KP_B855_130304_p1c2.mat'}

for iexpt = 1:size(repexpts,1)
    rootname{iexpt} = repexpts{iexpt}(1:19);
    
end
unq_expts = unique(rootname);

%% sort through expt directory and find expts
% that meet the given criteria
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


repexpts = [];
trials=5;
exptind=1;

for iexpt=1:size(exptnames,2)
    
    thisexpt=exptnames{iexpt};
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
                    expt.analysis.params.baselinewin(1)=basetimes(1)+basebegin;
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
            %
            %             [s,v] = GetVarDistrib(vmexpt,stimcond);
            %             sigdata{exptind} = s;
            %             vardata{exptind} = v;
            %             allstims{exptind} = stimcond;
            repexpts{exptind} = thisexpt;
            %             basewin{exptind} = expt.analysis.params.baselinewin;
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

%%

trials = 5;
signalDB = 80;
nresid_stim = [];
nresid_base = [];
skew_base = [];
skew_stim = [];
base_spk_vm = [];
stim_spk_vm = [];
stim_spk_thresh = [];
stim_spk_thresh = [];
stimind = 1;

PlotRinRsCm = 1;
plotSpikeShape = 0;
plotResidVm = 0;


%get RinRsCm for each actual unique expt
%would be worth it to go back to the raw data and re-do these
%because raw data has more data point which would help accuracy for Rs
%measurement and therefore Tau/Cm in particular
Rin = [];
Rs = [];
Rv = [];
TaoCell = [];
TaoV = [];
Cm = [];
expt_stepdata = [];
expt_baselineVm = [];
V_f = [];
Rv_dv = [];
unqind = 1;
exptind = 1;
exptunq_ind = [];
for iunq = 1:size(unq_expts,2)
    unq_expts{iunq};
    tmpstepdata = [];
    tmpbaselinedata = [];
    this_cell = unq_expts{iunq}
    for iexpt=1:size(repexpts,1)
        if ~isempty(regexp(repexpts{iexpt},unq_expts{iunq}))
            thisexpt=repexpts{iexpt};
            load([r.Dir.Expt thisexpt])
            stepdur = round(0.25/expt.wc.dt);
            stepstart = 553;
            vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
            tmpstepdata = [tmpstepdata;vmexpt.wc.data(:,stepstart:stepdur)*1000];
            tmpbaselinedata = [tmpbaselinedata; vmexpt.wc.data(:,1:stepstart)*1000];
            exptunq_ind (exptind) = iunq;
            exptind = exptind + 1;
        end
    end
    
    %only using trials that are in the lower 50% of the std values...
    %%%%%%%%%%%later should make sure that for each experiment there are
    %%%%%%%%%%%enough trials left at that point to get a good estimate...
    %%%%%%%%%%%%%%%%%%could just look at fits to see if reasonable?
    
    %%%%%%%**************
    %%%********DONT SELECT OUT LOWEST STD TRACES BECAUSE UNDERESTIMATING
    %%%HYPERPOL INDUCED CONDUCTANCE**************
    %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%
    %****NO TO THE ABOVE>>> TRIED THAT AND IT IS NOT NOT NOT BETTER
    std_distrib = std(tmpstepdata(:,300:round((stepdur-stepstart)/2))');
    %      figure;hist(std_distrib)
    keepinds = find(std_distrib <= median(std_distrib));
    tmpstepdata = tmpstepdata(keepinds,:);
    expt_stepdata(iunq,:) = mean(tmpstepdata);
    expt_baselineVm(iunq) = mean(mean(tmpbaselinedata));
    
    [out_struct, hfig,hfig2,hfig3] = MetaResponseAnal_RsRin(tmpstepdata,expt);
    allfields = fieldnames(out_struct);
    for ifield = 1:size(allfields,1)
        s = [allfields{ifield} '(unqind) = out_struct.' allfields{ifield} ';'];
        eval(s)
    end
    
    if PlotRinRsCm == 1
        foldername = '/Users/kperks/GitHub/Data_Mat/Analysis/Meta_Responses_DB/RinRsCm/';
        saveas(hfig,[foldername unq_expts{iunq} 'Stepdata.fig'])
        saveas(hfig,[foldername unq_expts{iunq} 'Stepdata.png'])
        
        if ~isempty(hfig2)
            saveas(hfig2,[foldername unq_expts{iunq} 'Stepdata_3sumfit.fig'])
            saveas(hfig2,[foldername unq_expts{iunq} 'Stepdata_3sumfit.png'])
        end
    end
    
    close(hfig)
    close(hfig2)
    close(hfig3)
    unqind = unqind + 1;
end

for iunq = 1:size(expt_stepdata,1)
    allsteps(iunq,:) = expt_stepdata(iunq,:) - expt_stepdata(iunq,1);
end
figure
line([1:size(allsteps,2)]*expt.wc.dt,allsteps,'color',[0.5 0.5 0.5])

for iunq = 1:size(expt_stepdata,1)
    allsteps(iunq,:) = expt_stepdata(iunq,:) - min(expt_stepdata(iunq,:));
    allsteps(iunq,:) = allsteps(iunq,:) / max(allsteps(iunq,:));
end
figure
line([1:size(allsteps,2)]*expt.wc.dt,allsteps,'color',[0.5 0.5 0.5],'LineWidth',1)
line([1:size(allsteps,2)]*expt.wc.dt,mean(allsteps),'color','k','LineWidth',3)

if PlotPop == 1
    edges = min(TaoCell):(max(TaoCell)-min(TaoCell))/5:max(TaoCell);
    bins = histc(TaoCell, edges);
    hfig = figure;
    bar(edges*1000,bins/20)
    xlabel('msec','FontSize',14)
    ylabel('proportion of cells with each Tao','FontSize',14)
    title('Time Constant n = 20','FontSize',22)
    
    edges = min(TaoV):(max(TaoV)-min(TaoV))/5:max(TaoV);
    bins = histc(TaoV, edges);
    hfig = figure;
    bar(edges*1000,bins/20)
    xlabel('msec','FontSize',14)
    ylabel('proportion of cells with each voltage-dependent Tao','FontSize',14)
    title('n = 20','FontSize',22)
    
    edges = min(Rin):(max(Rin)-min(Rin))/5:max(Rin);
    bins = histc(Rin, edges);
    hfig = figure;
    bar(edges,bins/20)
    xlabel('megaOhm','FontSize',14)
    ylabel('proportion of cells with each Rin','FontSize',14)
    title('Rin n = 20','FontSize',22)
    
    edges = min(Rs):(max(Rs)-min(Rs))/5:max(Rs);
    bins = histc(Rs, edges);
    hfig = figure;
    bar(edges,bins/20)
    xlabel('megaOhm','FontSize',14)
    ylabel('proportion of cells with each Rs','FontSize',14)
    title('Rs n = 20','FontSize',22)
    
    figure
    scatter([1:20], Rs)
    set(gca,'XTick',[1:20])
    set(gca,'XTickLabel',unq_expts)
    rotateXLabels(gca,90)
    line([1,20],[60,60],'color','k')
    line([1,20],[50,50],'color','k')
    set(gca,'XLim',[1,20])
    ylabel('access resistance','FontSize',14)
    
    figure;
    scatter(Rin,TaoCell*1000,100,'k','fill')
    ylabel('Time Constant (msec)','FontSize',14)
    xlabel('Input Resistance (megaOhm)','FontSize',14)
    
    figure;
    scatter(Rin,Cm,100,'k','fill')
    ylabel('Cell Capacitance (pFarad)','FontSize',14)
    xlabel('Input Resistance (megaOhm)','FontSize',14)
    
    figure;
    scatter(Rin,TaoV*1000,100,'k','fill')
    ylabel('Voltage-Dependent Time Constant (msec)','FontSize',14)
    xlabel('Input Resistance (megaOhm)','FontSize',14)
    
    figure;
    scatter(Rin,(Rv_dv),100,'k','fill')
    xlabel('Input Resistance','FontSize',14)
    ylabel('dv due to v-dependent g','FontSize',14)
    %this should actually depend on the value of the membrane potential at
    %the beginning of the voltage-dependent conductance
    
    figure;
    scatter(V_f,(Rv_dv),100,'k','fill')
    xlabel('steady-state v during step','FontSize',14)
    ylabel('dv due to v-dependent g','FontSize',14)
    
    
end
%
%     unqind = unqind + 1;
%    end
%
%
% end

%spike shapes
%don't care about requiring a particular dB or a particular number of reps

%%%%get Rs for each cell:signal group so that can filter for only low Rs
spkind = 1;
RsSpk = [];
for iexptverb = 1:size(exptunq_ind,2)
    unqind = exptunq_ind(iexptverb);
    RsVerbose (iexptverb) = Rs(unqind);
    VrestVerbose (iexptverb) = expt_baselineVm(unqind);
    if RsVerbose (iexptverb) < 50
        if iexptverb ~=22
            RsSpk (spkind) = unqind;
            VrestSpk (spkind) = expt_baselineVm(unqind);
            spkind = spkind +1;
        end
    end
end

clipped_vm = []; %this is calculated within the for loop below
spk_vm = []; %this is the output of MetaResponseAnal_Spikes
VmBins = [];
SpkVmBins = [];
spk_shape = [];
spk_dv = [];
spk_thrsh = [];
useSpks = [];
spkwid = [];
spkheight = [];
spkdvdt = [];

%********since for some cells the stimulus lengths are different,
%********cannot do this analysis for each cell lumped... need to do it for
%********each cell/signal pair and then lump them after...
% for iunq = 1:size(unq_expts,2)
%     unq_expts{iunq};
%     %restrict this analysis to cells that have a low access resistance
%     %this will minimize lowpass filtering
%     if Rs(iunq) < 50
%     allstepdata = [];
%     this_cell = unq_expts{iunq}
unqind = 1;
exptname_spk = [];
for iexpt=1:size(repexpts,1)
    
    %         if ~isempty(regexp(repexpts{iexpt},unq_expts{iunq}))
   
    if RsVerbose(iexpt) < 50  && iexpt ~=22
         exptname_spk  = [exptname_spk;repexpts(iexpt)];
        thisexpt=repexpts{iexpt};
        load([r.Dir.Expt thisexpt])
        vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
        table=getClampTab(expt,{'clamp',0});
        thiscond=getsubstimcond(expt.stimcond,table.sigsplayed);
        [sigon,sigoff]=GetSigTimes(expt,thiscond,1);
        
        basetimes = vmexpt.analysis.params.baselinewin;
        %             [sigon,sigoff]=GetSigTimes(vmexpt,thiscond,1);
        
        sigexpt = filtesweeps(vmexpt,0,'wavnames',table.sigsplayed);
        sigdata = sigexpt.wc.data*1000;
        sigdata_filt = medfilt1(sigdata,200,[],2);
        highpassdata=HighpassGeneral(sigdata,1/expt.wc.dt);
        
        %             outcell=PlotAndAsk(highpassdata,'spk_thresh','negative');
        %             cellfun(@eval,outcell);
        %             %if spikes are not big enough... continue and dond't use this expt
        %             if spk_thresh < 10
        %                 continue
        %             end
        spk_thresh = 10;
        negative = 0;
        
        %set up input_struct with data to pass analysis functions
        input_struct.dt = expt.wc.dt;
        input_struct.spk_thresh = spk_thresh;
        input_struct.sigdata = sigdata(:,basetimes(1):end);
        input_struct.sigdata_filt = sigdata_filt(:,basetimes(1):end);
        input_struct.highpassdata = highpassdata(:,basetimes(1):end);
        
        % get spike shapes and spike threshold and actual spike peak times to get vm during spike
        
        %********all of a sudden getting errors...
        [out_struct, hfig ] = MetaResponseAnal_Spikes(input_struct);
        allfields = fieldnames(out_struct);
        for ifield = 1:size(allfields,1)
            s = [allfields{ifield} '{unqind} = out_struct.' allfields{ifield} ';'];
            eval(s)
        end
        if ~isempty(hfig)
            if plotSpikeShape == 1;
                ylims = get(gca,'YLim');
                title_str = [expt.name '    #spikes = ' num2str(size(spk_shape{unqind},1))];
                
                text(-0.025,ylims(2)-10,title_str,'Interpreter','none')
                saveas(hfig(1),[r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikeShapes/' ...
                    expt.name '_allCellSpk.png'])
                saveas(hfig(1),[r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikeShapes/' ...
                    expt.name '_allCellSpk.fig'])
            end
            
            
            close(hfig)
        end
        
        %     if isempty(spk_shape{unqind})
        %
        %     end
        
        if ~isempty(spk_shape{unqind})
            
            input_struct.shape = spk_shape{unqind};
            input_struct.initVm = spk_initVm{unqind};
            
            [out_struct,hfig] = MetaResponseAnal_SpikeCluster(input_struct);
            
            allfields = fieldnames(out_struct);
            for ifield = 1:size(allfields,1)
                s = [allfields{ifield} '{unqind} = out_struct.' allfields{ifield} ';'];
                eval(s)
            end
            
            nextfig = size(hfig,2)+1;
            hfig(nextfig) = figure;
            hold on
            subplot(3,1,1)
            hist(spkwid{unqind})
            title('spike width (msec)')
            subplot(3,1,2)
            hist(spkheight{unqind})
            title('spike height (mV)')
            subplot(3,1,3)
            hist(spkdvdt{unqind})
            title('dvdt initial rise (V/sec)')
            set(hfig(nextfig),'Position',[98   131   463   663])
            
            % save figures
            if plotSpikeShape == 1;
                
                saveas(hfig(1),[r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikeShapes/' ...
                    expt.name '_SpkClusters.png'])
                saveas(hfig(1),[r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikeShapes/' ...
                    expt.name '_SpkClusters.fig'])
                
                saveas(hfig(2),[r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikeShapes/' ...
                    expt.name '_SpksFirstCluster.png'])
                saveas(hfig(2),[r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikeShapes/' ...
                    expt.name '_SpksFirstCluster.fig'])
                
                saveas(hfig(3),[r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikeShapes/' ...
                    expt.name '_SpksSecondCluster.png'])
                saveas(hfig(3),[r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikeShapes/' ...
                    expt.name '_SpksSecondCluster.fig'])
                
                saveas(hfig(4),[r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikeShapes/' ...
                    expt.name '_SpksDistAll.png'])
                saveas(hfig(4),[r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikeShapes/' ...
                    expt.name '_SpksDistAll.fig'])
                
                saveas(hfig(5),[r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikeShapes/' ...
                    expt.name '_SpksDistFirstClus.png'])
                saveas(hfig(5),[r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikeShapes/' ...
                    expt.name '_SpksDistFirstClus.fig'])
            end
            
            
            % close figures
            for ifig = 1:size(hfig,2)
                close(hfig(ifig))
            end
            %get "p(spike|Vm)" (sort of) for all spikes
            %do whole distribution of Vm
            %overlap that with distribution of Vm "under" each spike
        end
        VmEdges_spk = [-90:1:-10];
        input_struct.VmEdges = VmEdges_spk;
        input_struct.highpassdata = highpassdata(:,basetimes(1):end);
        input_struct.sigdata_filt = sigdata_filt(:,basetimes(1):end);
        out_struct = MetaResponseAnal_SpkVmDist(expt,input_struct);
        
        allfields = fieldnames(out_struct);
        for ifield = 1:size(allfields,1)
            s = [allfields{ifield} ' = out_struct.' allfields{ifield} ';'];
            eval(s)
        end
        all_VmBins{unqind} = VmBins;
        all_SpkVmBins{unqind} = SpkVmBins;
        all_numSpks{unqind} = numSpks;
        all_numVm{unqind} = numVm;
        all_pSpk{unqind} = pSpk;
        
        input_struct.highpassdata = highpassdata(:,basetimes(1):basetimes(2));
        input_struct.sigdata_filt = sigdata_filt(:,basetimes(1):basetimes(2));
        out_struct = MetaResponseAnal_SpkVmDist(expt,input_struct);
        
        allfields = fieldnames(out_struct);
        for ifield = 1:size(allfields,1)
            s = [allfields{ifield} ' = out_struct.' allfields{ifield} ';'];
            eval(s)
        end
        
        spont_VmBins{unqind} = VmBins;
        spont_SpkVmBins{unqind} = SpkVmBins;
        spont_numSpks{unqind} = numSpks;
        spont_numVm{unqind} = numVm;
        spont_pSpk{unqind} = pSpk;
        
        input_struct.highpassdata = highpassdata(:,sigon:sigoff);
        input_struct.sigdata_filt = sigdata_filt(:,sigon:sigoff);
        out_struct = MetaResponseAnal_SpkVmDist(expt,input_struct);
        
        allfields = fieldnames(out_struct);
        for ifield = 1:size(allfields,1)
            s = [allfields{ifield} ' = out_struct.' allfields{ifield} ';'];
            eval(s)
        end
        stim_VmBins{unqind} = VmBins;
        stim_SpkVmBins{unqind} = SpkVmBins;
        stim_numSpks{unqind} = numSpks;
        stim_numVm{unqind} = numVm;
        stim_pSpk{unqind} = pSpk;
        %in a different analysis can sort what spikes are spont and
        %what spikes are driven
        %might need specific data targeted at this where I get more
        %silence recordings
        
        unqind = unqind + 1;
    end
end
%     end
%     end

%get mean spike shape and params for each unique cell averaged over all
%stimuli (and experiements)
unq_spk_shape_norm = [];
unq_spk_shape = [];
unq_wid = [];
unq_height = [];
unq_dvdt = [];
unq_initVm = [];
unq_Vrest = [];
ii = 1;
for iunq = 1:size(unq_expts,2)
    unqinds = find(RsSpk == iunq);
    if ~isempty(unqinds)
        thisSpks = useSpks(unqinds);
        thiswid = spkwid(unqinds);
        thisheight = spkheight(unqinds);
        thisdvdt = spkdvdt(unqinds);
        thisinitVm = spk_initVm(unqinds);
        thisVrest = VrestSpk(unqinds);
        tmpshape = [];
        tmpwid = [];
        tmpheight = [];
        tmpdvdt = [];
        tmpinitVm = [];
        tmpVrest = [];
        for icond = 1:size(thisSpks,2)
            tmpshape =  [tmpshape; thisSpks{icond}];
            tmpwid =  [tmpwid, thiswid{icond}];
            tmpheight =  [tmpheight, thisheight{icond}];
            tmpdvdt =  [tmpdvdt, thisdvdt{icond}];
            tmpinitVm =  [tmpinitVm, thisinitVm{icond}];
            tmpVrest =  [tmpVrest, thisVrest(icond)];
        end
        
        unq_spk_shape(ii,:) = mean(tmpshape,1);
        unq_spk_shape_norm(ii,:) = mean(tmpshape,1) / max(mean(tmpshape,1));
        unq_wid(ii,:) = mean(tmpwid);
        unq_height(ii,:) = mean(tmpheight);
        unq_dvdt(ii,:) = mean(tmpdvdt);
        unq_initVm(ii,:) = mean(tmpinitVm);
        unq_initVmAll{ii} = tmpinitVm;
        unq_Vrest(ii,:) = mean(tmpVrest);
        ii = ii+1;
    end
end

figure;
line([1:size(unq_spk_shape_norm,2)]*expt.wc.dt*1000,unq_spk_shape_norm,'color',[0.5 0.5 0.5])
line([1:size(unq_spk_shape_norm,2)]*expt.wc.dt*1000,mean(unq_spk_shape_norm),...
    'color','k','LineWidth',5)
axis tight
set(gca,'YLim',[-0.3,1.1])
xlabel('msec','FontSize',14)
ylabel('normalized height','FontSize',14)
title(['n cells = ' num2str(size(unq_spk_shape,1))],'FontSize',14)


figure;
line([1:size(unq_spk_shape,2)]*expt.wc.dt*1000,unq_spk_shape,'color',[0.5 0.5 0.5])
line([1:size(unq_spk_shape,2)]*expt.wc.dt*1000,mean(unq_spk_shape),...
    'color','k','LineWidth',5)
axis tight
%set(gca,'YLim',[-0.3,1.1])
xlabel('msec','FontSize',14)
ylabel('height','FontSize',14)
title(['n cells = ' num2str(size(unq_spk_shape,1))],'FontSize',14)

hfig = figure;
hold on
subplot(3,1,1)
hist(unq_wid)
title('spike width (msec)')
subplot(3,1,2)
hist(unq_height)
title('spike height (mV)')
subplot(3,1,3)
hist(unq_dvdt)
title('dvdt initial rise (V/sec)')
set(hfig,'Position',[98   131   463   663])


X = [unq_wid,unq_height,unq_dvdt];
[idx,ctrs] = kmeans(X,2);

hfig= figure
scatter3(X(idx==1,1),X(idx==1,2),X(idx==1,3),12,'r','fill')
hold on
scatter3(X(idx==2,1),X(idx==2,2),X(idx==2,3),12,'b','fill')
plot3(ctrs(:,1),ctrs(:,2),ctrs(:,3),'kx',...
    'MarkerSize',12,'LineWidth',2)
plot3(ctrs(:,1),ctrs(:,2),ctrs(:,3),'ko',...
    'MarkerSize',12,'LineWidth',2)
xlabel('spkwid msec')
ylabel('spkheight')
zlabel('spkdvdt V/sec')
title(['n = ' num2str(size(idx,1))]);



%need a standard for spike height, overshoot, and variance of height
%can filter some based on Rs because i know that causes lowpass filtering
%i can only use cells that don't pass for measure of variance in Vm
ncols = 4;
[colvec colsize rowvec rowsize] = subplotinds(ncols,size(all_SpkVmBins,2));
VmEdges = [-90:1:-10];
hfig = figure;
hold on
for icell = 1:size(all_SpkVmBins,2)
    subplot('Position',[colvec(icell),rowvec(icell),colsize,rowsize])
    set(gca,'YLim',[0,1], 'XLim',[-60,-20])
    hold on
    stairs(VmEdges,spont_SpkVmBins{icell},'color','b','LineWidth',3)
    stairs(VmEdges,stim_SpkVmBins{icell},'color','r','LineWidth',3)
    %     stairs(VmEdges,spont_VmBins{icell},'color','b','LineWidth',3)
    %     stairs(VmEdges,stim_VmBins{icell},'color','r','LineWidth',3)
    %
    
    if icell==size(colvec,1)+1
        title('Vm distribution underlying spikes;   blue : pre-stim   red : Vm under spike','FontSize',14)
    end
    if icell ~= size(colvec,1)
        set(gca,'XTick',[],'YTick',[])
    end
    title(repexpts{icell},'Interpreter','none')
end
set(hfig,'Position',pos);



hfig = figure;
hold on
for icell = 1:size(all_SpkVmBins,2)
    
    subplot('Position',[colvec(icell),rowvec(icell),colsize,rowsize])
    set(gca,'XLim',[-80,-20])
    hold on
%     stairs(VmEdges,spont_pSpk{icell},'color','b','LineWidth',3)
%     stairs(VmEdges,stim_pSpk{icell},'color','r','LineWidth',3)
    stairs(VmEdges,stim_pSpk{icell}+spont_pSpk{icell},'color','k','LineWidth',3)
    
    if icell==size(colvec,1)+1
        title('p(spike) at each Vm;   black : allVm   red : Vm under spike','FontSize',14)
    end
    if icell == size(colvec,1)
        legend('pre-stimulus','during stimulus')
    end
    if icell ~= size(colvec,1)
        set(gca,'XTick',[])%,'YTick',[])
    end
    title(['bird: ' repexpts{icell}(5:7) '    cell: ' repexpts{icell}(16:19)])
end

set(hfig,'Position',pos);


hfig = figure;
hold on
for icell = 1:size(unq_initVmAll,2)
    
%     subplot('Position',[colvec(icell),rowvec(icell),colsize,rowsize])
%     set(gca,'XLim',[-80,-20])
    hold on
%     stairs(VmEdges,spont_pSpk{icell},'color','b','LineWidth',3)
%     stairs(VmEdges,stim_pSpk{icell},'color','r','LineWidth',3)
%     stairs(VmEdges,stim_pSpk{icell}+spont_pSpk{icell},'color','k','LineWidth',3)
[x,p] = empcdf(unq_initVmAll{icell});
     stairs(x,p,'color','k','LineWidth',3)
%     if icell==size(colvec,1)+1
%         title('p(spike) at each Vm;   black : allVm   red : Vm under spike','FontSize',14)
%     end
%     if icell == size(colvec,1)
%         legend('pre-stimulus','during stimulus')
%     end
%     if icell ~= size(colvec,1)
%         set(gca,'XTick',[])%,'YTick',[])
%     end
%     title(['bird: ' repexpts{icell}(5:7) '    cell: ' repexpts{icell}(16:19)])
end

set(hfig,'Position',pos);

%plot normalized spike shape for cells with low Rs
% and consistent spike shape
hfig = figure;
hold on
for icell = 1:size(all_SpkVmBins,2)
    normspk = [];
end

%% stims from Jason's data
% from batchbatchGetStimsionto.m
% allwavs_ionto
% allnames_ionto
a = cellstr(allnames_ionto);
[c,ia,ic_respmax]= unique(a);

for iwav = 1:size(allwavs,2)
    wavlen_ionto(iwav) = size(allwavs_ionto{iwav},1);
end

minlen_ionto = min(wavlen_ionto);
for iwav = 1:size(allwavs_ionto,2)
    clipwavs_ionto(iwav,:) = allwavs_ionto{iwav}(1:minlen_ionto,1)';
end

%actually don't want to get unique ones... because want a representation of
%power actually played, not power of each stim avgd once

% fs = FS;
% unqwaves_ionto = clipwavs_ionto(ia,:);

rmswav_ionto = [];
for iwav = 1:size(clipwavs_ionto,1)
    y = clipwavs_ionto(iwav,:);
    dcoff = (mean(y));
    nodc = y-dcoff;
    rmswav_ionto(iwav,:) = sqrt((nodc.^2));
end

meanfiltRMS_ionto = medfilt1(mean(rmswav_ionto),200,[],2);
normRMS_ionto = meanfiltRMS_ionto/max(meanfiltRMS_ionto);
figure;
hold on
line([1:size(normRMS_ionto,2)]/FS ,normRMS_ionto,'color','k','LineWidth',2)
normspk_base_clip = normspk_base(1,50:end-50) - median(normspk_base(1,50:2000));
line(([1:size(normspk_base_clip,2)]/10000)-2,normspk_base_clip,'color','r','LineWidth',2)
fanonorm = Result_base.FanoFactor - median(Result_base.FanoFactor(300:end,1));
fanonorm = fanonorm / max(fanonorm);
line(Result_base.times/1000,fanonorm,'color','b','LineWidth',2)
set(gca,'XLim',[-2,5])

%% get spikes mean psth across Jason's data


%run first section of batchbatchPSTHionto.m to get base_data_All
for icond = 1:max(size(base_data_All))
    tmpspks = zeros(size(base_data_All(icond).spikes,1),size(base_data_All(icond).spikes,2));
    for itrial = 1:size(base_data_All(icond).spikes,1)
        tmpspks(itrial,find(base_data_All(icond).spikes(itrial,:))) = 1;
    end
    trial_len(icond) = size(base_data_All(icond).spikes,2);
    meanspks{icond} = mean(tmpspks);
end

for icond = 1:max(size(base_data_All))
    spkvec_base(icond,:) = meanspks{icond}(1,1:min(trial_len));
end

gausstosmooth = fspecial('gaussian', [20,1],2);
gausstosmooth = gausstosmooth/max(gausstosmooth);
psth_base = conv(gausstosmooth,mean(spkvec_base));
psth_base = medfilt1(psth_base,200,[],2);
all_spk_base = [];
%get all spikes trials psth and smoothed to randomly choose trials to
%average and get cross correlation between
for icond = 1:size(spkvec_base,1)
    all_spk_base (icond,:) = conv(gausstosmooth,spkvec_base(icond,:));
end
filt_allspk_nomean = medfilt1(all_spk_base,200,[],2);
filt_allspk_meanfirst = medfilt1(mean(all_spk_base),200,[],2);


for icond = 1:max(size(gz_data_All))
    tmpspks = zeros(size(gz_data_All(icond).spikes,1),size(gz_data_All(icond).spikes,2));
    for itrial = 1:size(gz_data_All(icond).spikes,1)
        tmpspks(itrial,find(gz_data_All(icond).spikes(itrial,:))) = 1;
    end
    trial_len(icond) = size(gz_data_All(icond).spikes,2);
    meanspks{icond} = mean(tmpspks);
end

for icond = 1:max(size(gz_data_All))
    spkvec_gz(icond,:) = meanspks{icond}(1,1:min(trial_len));
end

gausstosmooth = fspecial('gaussian', [20,1],2);
gausstosmooth = gausstosmooth/max(gausstosmooth);
psth_gz = conv(gausstosmooth,mean(spkvec_gz));
psth_gz = medfilt1(psth_gz,200,[],2);


%%
%%%%% prep for Vm variance for each cell: signal pair
for iexpt=1:size(repexpts,1)
    %         if ~isempty(regexp(repexpts{iexpt},unq_expts{iunq}))
    thisexpt=repexpts{iexpt};
    load([r.Dir.Expt thisexpt])
    vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
    table=getClampTab(expt,{'clamp',0});
    thiscond=getsubstimcond(expt.stimcond,table.sigsplayed);
    [sigon,sigoff]=GetSigTimes(expt,thiscond,1);
    allon (iexpt) = sigon;
    alloff(iexpt) = sigoff;
    alldatalen(iexpt) = size(vmexpt.wc.data,2);
end

%% get stimulus wavs
sig_wavname = [];
sig_wav = [];
tmpnames = [];
sigind = 1;
for iexpt=1:size(repexpts,1)
    %         if ~isempty(regexp(repexpts{iexpt},unq_expts{iunq}))
    thisexpt=repexpts{iexpt};
    load([r.Dir.Expt thisexpt])
    vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
    table=getClampTab(expt,{'clamp',0});
    keepsigs=reprequire(table,trials);
    thiscond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
    
    for isig = 1:size(thiscond,2)
        sig_wav{sigind} = thiscond(isig).wavs;
        %         A = {'one','two','twenty-two','One','two'};
        siglen(sigind) = size(thiscond(isig).wavs,1)
        tmpnames = strvcat(tmpnames ,thiscond(isig).wavnames);
        sigind = sigind + 1;
    end
    %%%need to do this where I do an average for each cell...
    % instead of this where every cell:signal pair is its own average
end

a = cellstr(tmpnames);
[c,ia,ic_respmax]= unique(a);

allwavs = [];
minlen = min(siglen);
for iwav = 1:size(sig_wav,2)
    allwavs(iwav,:) = sig_wav{iwav}(1:minlen,1)';
end

fs = 44100;
unqwaves = allwavs(ia,:);
rmswav = [];
filtwav = [];
for iwav = 1:size(unqwaves,1)
    filtwav(iwav,:) = fftFilter(unqwaves(iwav,:)',fs,50,2)';
    y = filtwav(iwav,:);
    dcoff = (mean(y));
    nodc = y-dcoff;
    rmswav(iwav,:) = sqrt((nodc.^2));
end

meanRMSstims = medfilt1(mean(rmswav),200,[],2);
figure;
line([1:size(meanRMSstims,2)]/fs,meanRMSstims/max(meanRMSstims),'color','k','LineWidth',4)
normVm = mean(vm_mean) - min(mean(vm_mean));
normVm = normVm / max(normVm);
vmtime = [1:size(normVm,2)]*expt.wc.dt - 0.5;
line(vmtime,normVm,'color','r','LineWidth',4)
normVar = mean(vm_var) - min(mean(vm_var));
normVar = normVar / max(normVar);
line(vmtime,normVar,'color','b','LineWidth',4)


normvar = normVar(1,round(0.5/expt.wc.dt):end);
normvar = normvar/max(normvar);
olddt=1/fs;
bin=round(1/olddt/round(1/expt.wc.dt));
newdt=olddt*bin;
dnsampRMS=meanRMSstims(:,1:bin:end);
dnsampRMS = dnsampRMS/max(dnsampRMS);
minlen = min([size(dnsampRMS,2),size(normvar,2)]);
normvar = normvar(1,1:minlen);
dnsampRMS = dnsampRMS(1,1:minlen);
figure;
hold on
scatter(dnsampRMS,normvar,20,'k','fill')
cf = fit(dnsampRMS',normvar','exp1');
% fitline = cf.p1*dnsampRMS + cf.p2;
fitline = cf.a*exp(cf.b*dnsampRMS)
line(dnsampRMS,fitline,'color','r')
[r rmse] = rsquare(dnsampRMS,fitline);
xlabel('normalized RMS','FontSize',14)
ylabel('normalized Vm Variance','FontSize',14)
text(0.5,0.9,[num2str(cf.a) '*exp(' num2str(cf.b) '*RMS)'],'FontSize',14)

edges = [0:0.1:1];
[nvar,bin] = histc(dnsampRMS,edges);
for ibin = 1:size(nvar,2)
    bininds = find(bin == ibin);
    binRMS (ibin) = mean(dnsampRMS(bininds));
    binVAR (ibin) = mean(normvar(bininds));
end
figure;
hold on
scatter(binRMS,binVAR,100,'k','fill')
xlabel('mean RMS per bin','FontSize',14)
ylabel('mean variance per binned mean RMS','FontSize',14)
set(gca,'YLim',[0,1],'XLim',[0,1])
cf = fit(binRMS',(binVAR - min(binVAR))','exp1');
fitline = cf.a*exp(cf.b*binRMS) + min(binVAR);
line(binRMS,fitline,'color','r')
[r rmse] = rsquare(binVAR,fitline);
title(['r-square for exp1 fit = ' num2str(r)...
    '  ;  offset (min variance) = ' num2str(min(binVAR))],'FontSize',18)

figure;
hold on
line([1:size(dnsampRMS,2)]*newdt,dnsampRMS,'color','r')
line([1:size(meanRMSstims,2)]/44100,meanRMSstims,'color','k')

%time course of stimulus rms increase and variance decrease
%go back to using "original" rms calc instead of medfilt to avoid artifacts
%of filtering using median at edges
options = fitoptions('exp2');

%fit stimulus RMS mean
endlen = round(size(rmswav,2)/4);
endt = endlen / fs; %seconds
meanRMS = mean(rmswav);
meanRMS = meanRMS / max(meanRMS);
endl = round(endt*fs);
xtime = [1:endl]/fs;
yRMS = meanRMS(1,1:endl);
tmp = yRMS - median(yRMS);
options.Lower = [-Inf,-Inf];
options.Upper = [0,0];
f = fit(xtime', tmp' , 'exp1' ,options);
fitline = f.a*exp(f.b*xtime);
fitline = fitline - min(fitline);
taoRMS = -1/f.b;
Tsamp = round(taoRMS * 44100);
[rD rmse] = rsquare(yRMS(1,1:Tsamp),fitline(1,1:Tsamp));
figure;
line(xtime,yRMS,'color','k');
line(xtime,fitline,'color','r');

%fit Vm variance mean
endl = round(endt/expt.wc.dt);
xtime = [1:endl]*expt.wc.dt;
yVAR = normvar(1,1:end);
figure;
plot(yVAR');
xonset = input('xtime onset')
yVAR = normvar(1,xonset:endl+xonset-1);
tmp = yVAR - median(yVAR);
options.Lower = [0,-Inf];
options.Upper = [Inf,0];
f = fit(xtime', tmp' , 'exp1' ,options);
fitline = f.a*exp(f.b*xtime);
fitline = fitline + median(yVAR);
taoVAR = -1/f.b;
Tsamp = round(taoVAR /expt.wc.dt);
[rD rmse] = rsquare(yVAR(1,1:Tsamp),fitline(1,1:Tsamp));
figure;
line(xtime,yVAR,'color','k');
line(xtime,fitline,'color','r');

%fit Vm mean
normean = normVm(1,round(0.5/expt.wc.dt):end);
normean = normean/max(normean);
normean = normean(1,1:minlen);
endl = round(endt/expt.wc.dt);
xtime = [1:endl]*expt.wc.dt;
% yMEAN = normean(1,1:end);
figure;
plot(normVm');
xonset = input('xtime onset')
xonset = xonset - round(0.5/expt.wc.dt);
yMEAN = normean(1,xonset:endl+xonset-1);
tmp = yMEAN - median(normean(1,2000:7000));
options.Lower = [-Inf,-Inf];
options.Upper = [0,0];
f = fit(xtime', tmp', 'exp1' ,options);
fitline = f.a*exp(f.b*xtime);
fitline = fitline + median(normean(1,2000:7000));
taoMEAN = -1/f.b;
Tsamp = round(taoMEAN /expt.wc.dt);
[rD rmse] = rsquare(yMEAN(1,1:Tsamp),fitline(1,1:Tsamp));
figure;
line(xtime,yMEAN,'color','k');
line(xtime,fitline,'color','r');

VarDecParam.taoVAR = taoVAR;
VarDecParam.taoRMS = taoRMS;
VarDecParam.taoMEAN = taoVAR;
VarDecParam.VarDecOnset = xonset * expt.wc.dt;



%% from stimulus onset
vm_var = [];
vm_mean = [];
siglen = min(alloff - allon);
sigind = 1;
prestim_time = 1; %seconds
prestim_len = round(prestim_time/expt.wc.dt);
for iexpt=1:size(repexpts,1)
    %         if ~isempty(regexp(repexpts{iexpt},unq_expts{iunq}))
    thisexpt=repexpts{iexpt};
    load([r.Dir.Expt thisexpt])
    vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
    table=getClampTab(expt,{'clamp',0});
    keepsigs=reprequire(table,trials);
    thiscond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
    startind = allon(iexpt) - prestim_len;
    stopind = allon(iexpt) + siglen;
    
    for isig = 1:size(thiscond,2)
        sigexpt = filtesweeps(vmexpt,0,'wavnames',thiscond(isig).wavnames);
        sigdata = sigexpt.wc.data*1000;
        sigdata_filt = medfilt1(sigdata,200,[],2);
        
        vm_var(sigind,:) = var(sigdata_filt(:,startind:stopind));
        vm_mean(sigind,:) = mean(sigdata_filt(:,startind:stopind));
        sigind = sigind + 1;
    end
    %%%need to do this where I do an average for each cell...
    % instead of this where every cell:signal pair is its own average
end

%get ci on var
confint = [];
for isamp = 1:size(vm_var,2)
    confint(isamp,:) = getCDFconf (vm_var(:,isamp),95);
end

xtime = ([1:size(vm_var,2)]*expt.wc.dt)-prestim_time;
figure;
line(xtime,mean(vm_var),'color','k','LineWidth',3)
line(xtime, confint(:,1),'color',[0.5 0.5 0.5])
line(xtime, confint(:,2),'color',[0.5 0.5 0.5])
axis tight
set(gca,'YLim',[-10,90]);
ylims = get(gca,'YLim');
line([0,0],[ylims(1),ylims(2)],'LineStyle','--','color','k');
ylabel('mean variance ; 95% CI','FontSize',14)
xlabel('seconds','FontSize',14)

figure;
line(xtime,mean(vm_mean),'color','k','LineWidth',3)
axis tight
ylims = get(gca,'YLim');
line([0,0],[ylims(1),ylims(2)],'LineStyle','--','color','k');
ylabel('mean membrane potential (mV)','FontSize',14)
xlabel('seconds','FontSize',14)

%calculate when var drops below mean var during pre-stim...
%select that as start point to calculate time constant of decline of var

prestim_vars = mean(vm_var(:,1:(prestim_len)));
[xb, prestim_vardist] =empcdf(prestim_vars);
prestim_ci = getCDFconf (prestim_vars,99);

for isamp = (prestim_len):size(vm_var,2)
    thisvar = mean(vm_var(:,isamp));
    if thisvar < (prestim_ci(1))
        lowvar(isamp) = 1;
    else lowvar(isamp) = 0;
    end
end

winlen = 5;
for isamp = 1:size(lowvar,2) - winlen
    meanlen(isamp) = sum(lowvar(1,isamp:isamp+winlen));
end
signif_drop = crossing(meanlen,[],4);

% fitend = round(0.5/expt.wc.dt);
vmdrop = mean(vm_var(:,(min(signif_drop)):end));
% vmdrop = vmdrop(1,1:fitend) - median(vmdrop);
vmdrop = vmdrop - median(vmdrop);
vmdrop_time = [1:size(vmdrop,2)]*expt.wc.dt;
cf = fit(vmdrop_time',vmdrop','exp1');
fitline = cf.a*exp(cf.b*vmdrop_time);

figure;hold on
line(vmdrop_time,vmdrop,'color','k')
line(vmdrop_time,fitline,'color','r')

timeconstant = -1/cf.b


%%
% compare mean_vm with mean spike psth
spksamprate = 10000;

vmtime = ([1:size(vm_mean,2)]*expt.wc.dt)-prestim_time;
normvm = mean(vm_mean);
normvm = normvm - min(normvm(1,200:end-200));
normvm = normvm / abs(max(max(normvm)));

spktime = ([1:size(psth_base,2)]/spksamprate)-(2);
normspk_base = psth_base - min(psth_base(1,200:end-200));
normspk_base = normspk_base / max(normspk_base);

figure;
hold on
line(spktime, normspk_base','color','r');
line(vmtime, normvm','color','k');
axis tight
ylims = get(gca,'YLim');
line([0,0],[ylims(1),ylims(2)],'LineStyle','--','color','k')
legend('spikes from JVT data','Vm from my 20cells')
ylabel('normalized magnitude','FontSize',14)
xlabel('seconds','FontSize',14)

figure;imagesc(spktime,1,normspk_base)
colormap('gray')

% GZ condition: compare mean_vm with mean spike psth
spksamprate = 10000;

vmtime = ([1:size(vm_mean,2)]*expt.wc.dt)-prestim_time;
normvm = mean(vm_mean);
normvm = normvm - min(normvm);
normvm = normvm / abs(max(max(normvm)));

spktime = ([1:size(psth_gz,2)]/spksamprate)-(2);
normspk_gz = psth_gz - min(psth_gz(1,200:2000));
normspk_gz = normspk_gz / max(normspk_gz);

figure;
hold on
line(spktime, normspk_gz','color','r');
line(vmtime, normvm','color','k');
axis tight
ylims = get(gca,'YLim');
line([0,0],[ylims(1),ylims(2)],'LineStyle','--','color','k')
title('gabazine condition','FontSize',20)
legend('spikes from JVT data','Vm from my 20cells')
ylabel('normalized magnitude','FontSize',14)
xlabel('seconds','FontSize',14)


figure;
hold on
line(spktime, normspk_gz','color','r');
line(spktime, normspk_base'','color','k');
axis tight
ylims = get(gca,'YLim');
line([0,0],[ylims(1),ylims(2)],'LineStyle','--','color','k')
title('JVT data - mean spikes across population','FontSize',20)
legend('gabazine','baseline')
ylabel('normalized magnitude','FontSize',14)
xlabel('seconds','FontSize',14)

%
% % need to downsample
% spksamprate = 1000;
%
% vmtime = ([1:size(vm_var,2)]*expt.wc.dt)-prestim_time;
% dnsamp_xtime = dnsample_data(vmtime,1/expt.wc.dt,spksamprate);
% dnsamp_vm=dnsample_data(vm_mean,1/expt.wc.dt,spksamprate);
% normvm = dnsamp_vm / abs(max(max(dnsamp_vm)));
%
% spktime = ([1:size(psth_base,2)]/spksamprate)-(2000/spksamprate);
% normspk = psth_base / max(psth_base);
%
% figure;
% hold on
% line(dnsamp_xtime, (mean(normvm) - min(mean(normvm)))')
% line(spktime, (normspk - min(normspk(1,5:end-5)))')
%% from stimulus offset
off_vm_var = [];
off_vm_mean = [];
sigind = 1;
preoffset_len = min(alloff - allon);
postoffset_len = min(alldatalen - alloff);
for iexpt=1:size(repexpts,1)
    %         if ~isempty(regexp(repexpts{iexpt},unq_expts{iunq}))
    thisexpt=repexpts{iexpt};
    load([r.Dir.Expt thisexpt])
    vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
    table=getClampTab(expt,{'clamp',0});
    thiscond=getsubstimcond(expt.stimcond,table.sigsplayed);
    startind = alloff(iexpt)-preoffset_len;
    stopind = alloff(iexpt) + postoffset_len;
    
    for isig = 1:size(thiscond,2)
        sigexpt = filtesweeps(vmexpt,0,'wavnames',thiscond(isig).wavnames);
        sigdata = sigexpt.wc.data*1000;
        sigdata_filt = medfilt1(sigdata,200,[],2);
        
        off_vm_var(sigind,:) = var(sigdata_filt(:,startind:stopind));
        off_vm_mean(sigind,:) = mean(sigdata_filt(:,startind:stopind));
        sigind = sigind + 1;
    end
    %%%need to do this where I do an average for each cell...
    % instead of this where every cell:signal pair is its own average
end

%get ci on var
for isamp = 1:size(off_vm_var,2)
    confint(isamp,:) = getCDFconf (off_vm_var(:,isamp),95);
end

xtime = ([1:size(off_vm_var,2)]*expt.wc.dt)-preoffset_len*expt.wc.dt;
figure;
line(xtime,mean(off_vm_var),'color','k','LineWidth',3)
line(xtime, confint(:,1),'color',[0.5 0.5 0.5])
line(xtime, confint(:,2),'color',[0.5 0.5 0.5])
axis tight
set(gca,'YLim',[-10,90]);
ylims = get(gca,'YLim');
line([0,0],[ylims(1),ylims(2)],'LineStyle','--','color','k');
ylabel('mean variance ; 95% CI','FontSize',14)
xlabel('seconds','FontSize',14)

figure;
line(xtime,mean(off_vm_mean),'color','k','LineWidth',3)
axis tight
ylims = get(gca,'YLim');
line([0,0],[ylims(1),ylims(2)],'LineStyle','--','color','k');
ylabel('mean membrane potential (mV)','FontSize',14)
xlabel('seconds','FontSize',14)

%calculate when var drops below mean var during pre-stim...
%select that as start point to calculate time constant of decline of var

prestim_vars = mean(off_vm_var(:,1:(prestim_len)));
[xb, prestim_vardist] =empcdf(prestim_vars);
prestim_ci = getCDFconf (prestim_vars,99);

for isamp = (prestim_len):size(off_vm_var,2)
    thisvar = mean(off_vm_var(:,isamp));
    if thisvar < (prestim_ci(1))
        lowvar(isamp) = 1;
    else lowvar(isamp) = 0;
    end
end

winlen = 10;
for isamp = 1:size(lowvar,2) - winlen
    meanlen(isamp) = sum(lowvar(1,isamp:isamp+winlen));
end
signif_drop = crossing(meanlen,[],9);

% fitend = round(0.5/expt.wc.dt);
vmdrop = mean(off_vm_var(:,(min(signif_drop)):end));
% vmdrop = vmdrop(1,1:fitend) - median(vmdrop);
vmdrop = vmdrop - median(vmdrop);
vmdrop_time = [1:size(vmdrop,2)]*expt.wc.dt;
cf = fit(vmdrop_time',vmdrop','exp1');
fitline = cf.a*exp(cf.b*vmdrop_time);

figure;hold on
line(vmdrop_time,vmdrop,'color','k')
line(vmdrop_time,fitline,'color','r')


%%
%get Vm response windows and p(spike) in those windows?



VmRespVec_up = [];
VmRespVec_low = [];
VmBaseConf = [];
VmAllLen = [];
unqind = 1;
for iunq = 1:size(unq_expts,2)
    unq_expts{iunq};
    siglen = 0;
    allstepdata = [];
    this_cell = unq_expts{iunq};
    Vrest = [];
    response_vm_confint = [];
    response_vm_mean = [];
    response_vm_max = [];
    up_win = [];
    low_win = [];
    allVmRespVec = [];
    sigind = 1;
    for iexpt=1:size(repexpts,1)
        if ~isempty(regexp(repexpts{iexpt},this_cell))
            thisexpt=repexpts{iexpt};
            load([r.Dir.Expt thisexpt])
            
            vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
            % for now do not care if db or not
            %             [dbstimcond,dblevels]=getDBstimcond(vmexpt);
            table=getClampTab(expt,{'clamp',0});
            
            keepsigs=reprequire(table,trials);
            allsig=table.sigsplayed;
            repsigs=allsig(keepsigs);
            stimcond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
            
            for istim = 1:size(stimcond,2)
                sigexpt = filtesweeps(vmexpt,0,'wavnames',stimcond(istim).wavnames);
                input_struct.sigdata = sigexpt.wc.data*1000;
                input_struct.sigdata_filt = medfilt1(input_struct.sigdata,200,[],2);
                
                basetimes = sigexpt.analysis.params.baselinewin;
                [sigon,sigoff]=GetSigTimes(sigexpt,stimcond,istim);
                %set up input_struct with data to pass analysis functions
                input_struct.basetimes = basetimes;
                input_struct.sigon = sigon;
                input_struct.sigoff = sigoff;
                
                siglen = siglen + size(input_struct.sigdata_filt(:,basetimes(1):end),2);
                
                %prep spiking to be continuous and get spiking windows
                highpassdata=HighpassGeneral(input_struct.sigdata,1/expt.wc.dt); 
                [spikesmat, gausstosmooth]=getspikesmat(highpassdata,threshold,dt);
                
                
                % get vm response windows
                out_struct = MetaResponseAnal_VmResponseWin(expt,input_struct)
                allfields = fieldnames(out_struct);
                for ifield = 1:size(allfields,1)
                    s = [allfields{ifield} '{sigind} = out_struct.' allfields{ifield} ';'];
                    eval(s)
                end
                respfig = figure;
                hold on
                scaleticks = 1;
                ydatabound = [min(mean(input_struct.sigdata_filt(:,basetimes(1):end))), max(mean(input_struct.sigdata_filt(:,basetimes(1):end)))];
                xtime=[1:size(input_struct.sigdata_filt,2)]*expt.wc.dt;
                line(xtime,mean(input_struct.sigdata_filt),'color','k','LineWidth',3);
                plot([xtime(1),xtime(end)],[response_vm_confint{sigind}(2),response_vm_confint{sigind}(2)],'--','color','k')
                plot([xtime(1),xtime(end)],[response_vm_confint{sigind}(1),response_vm_confint{sigind}(1)],'--','color','k')
                SigTimeBox(gca, (sigon)*expt.wc.dt,sigoff*expt.wc.dt, get(gca,'YLim'),[0.5 0.5 0.5]);
                for iresp=1:size(up_win{sigind},2)
                    SigTimeBox(gca, up_win{sigind}(1,iresp)*expt.wc.dt, ...
                        up_win{sigind}(2,iresp)*expt.wc.dt, get(gca,'YLim'),'r');
                end
                for inhib=1:size(low_win{sigind},2)
                    SigTimeBox(gca, low_win{sigind}(1,inhib)*expt.wc.dt, ...
                        low_win{sigind}(2,inhib)*expt.wc.dt, get(gca,'YLim'),'b');
                end
                axis tight
                Vrest{sigind} = mean(min(input_struct.sigdata_filt(:,basetimes(1):basetimes(2))'));
                ylims = [Vrest{sigind},ydatabound(2)];
                xlims = get(gca,'XLim');
                xlims = [basetimes(1)*expt.wc.dt,xlims(2)];
                set(gca,'XLim',xlims);
                
                set(gca,'YLim',ylims,'YTick',...
                    [(floor(Vrest{sigind})-mod(floor(Vrest{sigind}),5)):5:(ceil(ymax)+mod(ceil(ymax),5))],...
                    'XTick',[0:1:floor(xtime(end))],'TickDir','out')
                text(xtime(1,basetimes(1)),round(Vrest{sigind}),[num2str(round(Vrest{sigind})) 'mV'],...
                    'HorizontalAlignment','center',	'BackgroundColor', 'k',...
                    'color',[1,1,1]);
                box off
                set(respfig,'Position',[212         523        1168         283])
                title([expt.name stimcond(istim).wavnames],'Interpreter','none')
                
                saveas(respfig, ...
                    ['/Users/kperks/GitHub/Data_Mat/Analysis/Meta_Responses_DB/VmRespWin/' ...
                    expt.name '_' stimcond(istim).wavnames '.fig']);
                saveas(respfig, ...
                    ['/Users/kperks/GitHub/Data_Mat/Analysis/Meta_Responses_DB/VmRespWin/' ...
                    expt.name '_' stimcond(istim).wavnames '.tif']);
                close(respfig)
                
                respfig = figure;
                hold on
                scaleticks = 1;
                ydatabound = [min(mean(input_struct.sigdata_filt(:,basetimes(1):end))), max(mean(input_struct.sigdata_filt(:,basetimes(1):end)))];
                xtime=[1:size(input_struct.sigdata_filt,2)]*expt.wc.dt;
                line(xtime,mean(input_struct.sigdata_filt),'color','k','LineWidth',3);
                plot([xtime(1),xtime(end)],[response_vm_confint{sigind}(2),response_vm_confint{sigind}(2)],'--','color','k')
                plot([xtime(1),xtime(end)],[response_vm_confint{sigind}(1),response_vm_confint{sigind}(1)],'--','color','k')
                line([ (sigon)*expt.wc.dt,sigoff*expt.wc.dt], [ydatabound(1)-1, ydatabound(1)-1],...
                    'color','k','LineWidth',3);
                for iresp=1:size(up_win{sigind},2)
                    line([up_win{sigind}(1,iresp)*expt.wc.dt, ...
                        up_win{sigind}(2,iresp)*expt.wc.dt],[ydatabound(1)-1.2, ydatabound(1)-1.2],...
                        'color' ,'r','LineWidth',0.5);
                end
                for inhib=1:size(low_win{sigind},2)
                    line([low_win{sigind}(1,inhib)*expt.wc.dt, ...
                        low_win{sigind}(2,inhib)*expt.wc.dt],[ydatabound(1)-1.4, ydatabound(1)-1.4],...
                        'color' ,'b','LineWidth',0.5);
                end
                
                axis tight
                Vrest{sigind} = mean(min(input_struct.sigdata_filt(:,basetimes(1):basetimes(2))'));
                ylims = [ydatabound(1)-2,ydatabound(2)];
                xlims = get(gca,'XLim');
                xlims = [basetimes(1)*expt.wc.dt,xlims(2)];
                set(gca,'XLim',xlims);
                
                set(gca,'YLim',ylims,'YTick',...
                    [(floor(Vrest{sigind})-mod(floor(Vrest{sigind}),5)):5:(ceil(ymax)+mod(ceil(ymax),5))],...
                    'XTick',[0:1:floor(xtime(end))],'TickDir','out')
                text(xtime(1,basetimes(1)),round(Vrest{sigind}),[num2str(round(Vrest{sigind})) 'mV'],...
                    'HorizontalAlignment','center',	'BackgroundColor', 'k',...
                    'color',[1,1,1]);
                box off
                                set(respfig,'Position',[212         523        1168         283])

                title([expt.name stimcond(istim).wavnames],'Interpreter','none')
                
                saveas(respfig, ...
                    ['/Users/kperks/GitHub/Data_Mat/Analysis/Meta_Responses_DB/VmRespWin/foreps_' ...
                    expt.name '_' stimcond(istim).wavnames '.fig']);
                saveas(respfig, ...
                    ['/Users/kperks/GitHub/Data_Mat/Analysis/Meta_Responses_DB/VmRespWin/foreps_' ...
                    expt.name '_' stimcond(istim).wavnames '.tif']);
                close(respfig)
                
                
                sigind = sigind + 1;
            end
        end
    end
    
    %summary data for each expt (cell/signal block)
    %VmVector of all Vm values during responses (up or down)
    tmpvm_up = [];
    tmpvm_low = [];
    tmpconfint = [];
    tmpvrest = [];
    
      for isig = 1:size(Vrest,2)
        tmpvrest = [tmpvrest; Vrest{isig}];
    end
        
    for isig = 1:size(allVmRespVec_up,2)
        tmpvm_up = [tmpvm_up, allVmRespVec_up{isig}];
    end
    for isig = 1:size(allVmRespVec_low,2)
          tmpvm_low = [tmpvm_low, allVmRespVec_low{isig}];
    end
    for isig = 1:size(response_vm_confint,2)
        tmpconfint = [tmpconfint; response_vm_confint{isig}];
    end
    VmRespVec_up{iunq} = tmpvm_up;
    VmRespVec_low{iunq} = tmpvm_low;
    VmBaseConf{iunq} = tmpconfint;
    VmAllLen(iunq) = siglen;
    VrestCell{iunq} = tmpvrest;
    
end

figure;
hold on
for iunq = 1:size(VmAllLen,2)
    [x,p]=empcdf(VmRespVec_up{iunq} - VmBaseConf(iunq,2));
    stairs(x,p,'color','r')
end


figure;
hold on
for iunq = 1:size(VmAllLen,2)
    [x,p]=empcdf(VmRespVec_up{iunq} - VmBaseConf(iunq,1));
      stairs(x,p,'color','b')
end

figure;
hold on
for iunq = 1:size(VmAllLen,2)
    [x,p]=empcdf(VmRespVec_up{iunq});
    stairs(x,p,'color','r')
    [x,p]=empcdf(VmRespVec_low{iunq});
     stairs(x,p,'color','b')
     
end
[x,p] = empcdf(unq_Vrest);
stairs(x,p,'color','k')
[x,p] = empcdf(unq_initVm);
stairs(x,p,'color','k')

n = [];
VmEdges = [-80:5:-20];
figure;
hold on
for iunq = 1:size(VmAllLen,2)
    n_vm_up(iunq,:) = histc(VmRespVec_up{iunq},VmEdges)/size(VmRespVec_up{iunq},2);
    n_vm_low(iunq,:) = histc(VmRespVec_low{iunq},VmEdges)/size(VmRespVec_low{iunq},2);
end
n_spk = histc(unq_initVm,VmEdges)/size(unq_initVm,1);
spkbins = find(n_spk);
mean_vmhist_up = mean(n_vm_up);
mean_vmhist_low = mean(n_vm_low);
n_rest = histc(unq_Vrest,VmEdges)/size(unq_Vrest,1);
stairs(VmEdges,n_rest,'color',[0.5 0.5 0.5],'LineWidth',4)
bar(VmEdges(1,spkbins)+2.5,mean_vmhist_up(1,spkbins),'y','BarWidth',1)
stairs(VmEdges,mean_vmhist_low,'color','b','LineWidth',4)
stairs(VmEdges,mean_vmhist_up,'color','r','LineWidth',4)
stairs(VmEdges,n_spk,'color','k','LineWidth',4)
title(['% of depol Vm response values overlapping with spike threshold range: ' ...
    num2str(round(sum(mean_vmhist_up(1,spkbins))*100))],'FontSize',14)


% get propoertion of Vm that is counted as part of
% depol response and proportion that is hyperpol response
for iunq = 1:size(VmAllLen,2)
    p_up(iunq) = size(VmRespVec_up{iunq},2)/VmAllLen(iunq);
    p_low(iunq) = size(VmRespVec_low{iunq},2)/VmAllLen(iunq);
end
figure;
hold on
scatter(repmat(1,size(p_low,1),size(p_low,2)),p_low - p_up)
% scatter(repmat(2,size(p_up,1),size(p_up,2)),p_up)


%%
    %     input_struct.spk_thresh = spk_thresh;
    
    
    %         %get power spectrum
    %         out_struct = MetaResponseAnal_PDS(expt, input_struct)
    %         allfields = fieldnames(out_struct);
    %         for ifield = 1:size(allfields,1)
    %             s = [allfields{ifield} ' = out_struct.' allfields{ifield} ';'];
    %             eval(s)
    %         end
    %
    %         %get residuals of Vm distributions for baseline and stimulus period
    %         MetaResponseAnal_VmResiduals(expt,input_struct)
    %         allfields = fieldnames(out_struct);
    %         for ifield = 1:size(allfields,1)
    %             s = [allfields{ifield} ' = out_struct.' allfields{ifield} ';'];
    %             eval(s)
    %         end
    %
    %
    VmEdges = [0:5:40];
    figure
    hold on
    for iexpt = 1:size(allVmRespVec,2)
        [n,bins] = histc(allVmRespVec{iexpt}-Vrest{iexpt},VmEdges);
        stairs(VmEdges,n/size(allVmRespVec{iexpt},2))
    end
    
    
    figure
    hold on
    for iexpt = 1:size(allVmRespVec,2)
        dv = allVmRespVec{iexpt}-Vrest{iexpt};
        
        [x,p] = empcdf(dv);
        %    [x,p] = empcdf(allVmRespVec{iexpt}-Vrest{iexpt});
        stairs(x,p)
    end
    
    figure
    hold on
    for iexpt = 1:size(allVmRespVec,2)
        
        
        [x,p] = empcdf(allVmRespVec{iexpt});
        %    [x,p] = empcdf(allVmRespVec{iexpt}-Vrest{iexpt});
        stairs(x,p)
    end
    
    figure
    hold on
    for iexpt = 1:size(allVmRespVec,2)
        dv = allVmRespVec{iexpt}-Vrest{iexpt};
        X = [dv;dv];
        [ic,ix] = kmeans(X',2);
        
        hfig= figure;
        scatter(dv(1,ic==1),dv(1,ic==1),12,'r','fill')
        hold on
        scatter(dv(1,ic==2),dv(1,ic==2),12,'b','fill')
        line([xax(1),xax(end)],[ix_respmax{sigind}(1),ix_respmax{sigind}(1)],'color','k','LineWidth',2)
        line([xax(1),xax(end)],[ix_respmax{sigind}(2),ix_respmax{sigind}(2)],'color','k','LineWidth',2)
        xlabel('response ind')
        ylabel('max mV')
        
        up_resp = dv(find(ic == (find(ix == max(ix)))));
        [x,p] = empcdf(up_resp);
        %    [x,p] = empcdf(allVmRespVec{iexpt}-Vrest{iexpt});
        stairs(x,p)
    end
    
    
    %%
    
    %%%%%%%%response windows and responses across DB
    % respfig = figure;
    % hold on
    % scaleticks = 1;
    % [grad,im]=colorGradient([0.08,0,0.5],[0,0.78,0.78],3);
    % ymax = max(max(ydatabound));
    % y_raster = ymax;
    % ylims(1) = min(min(ydatabound));
    % xtime=[1:size(vmexpt.wc.data,2)]*expt.wc.dt;
    % for idb=1:size(thisdb,2)
    % line(xtime(1,basetimes(1):end),sdata(idb,basetimes(1):end),'color',grad(idb,:),'LineWidth',3);
    % %             plot([xtime(1),xtime(end)],[confint(idb,2),confint(2)],'--','color','k')
    % %add a line between each raster plot
    % line([xtime(1),xtime(end)],[y_raster,y_raster],'color','k');
    % text(xtime(1),y_raster,[num2str(thisdb(idb)) 'dB SPL'],...
    % 'HorizontalAlignment','center',	'BackgroundColor', grad(idb,:),...
    % 'color',[1,1,1]);
    % for itrial=1:  nreps(idb)
    % for ispike=1:size(spiketimes{idb,itrial},2)
    % plot([spiketimes{idb,itrial}(ispike)*expt.wc.dt,...
    % spiketimes{idb,itrial}(ispike)*expt.wc.dt],...
    % [(itrial*scaleticks)+y_raster, ...
    % (scaleticks*(itrial+0.9))+y_raster], ...
    % 'color','k','LineWidth',2)
    % end
    % end
    % y_raster = y_raster + nreps (idb) +1;
    % end
    % axis tight
    % ylims = get(gca,'YLim');
    % SigTimeBox(gca, sigon*expt.wc.dt,sigoff*expt.wc.dt, get(gca,'YLim'),[0.5 0.5 0.5]);
    % for iresp=1:size(up_win,2)
    % SigTimeBox(gca, up_win(1,iresp)*expt.wc.dt, ...
    % up_win(2,iresp)*expt.wc.dt, get(gca,'YLim'),'r');
    % end
    % for inhib=1:size(low_win,2)
    % SigTimeBox(gca, low_win(1,inhib)*expt.wc.dt, ...
    % low_win(2,inhib)*expt.wc.dt, get(gca,'YLim'),'b');
    % end
    % set(gca,'YLim',[ylims(1),ylims(2)],'YTick',...
    % [(floor(ylims(1))-mod(floor(ylims(1)),5)):5:(ceil(ymax)+mod(ceil(ymax),5))],...
    % 'XTick',[0:1:floor(xtime(end))],'TickDir','out')
    % text(xtime(1),round(ylims(1)),[num2str(round(ylims(1))) 'mV'],...
    % 'HorizontalAlignment','center',	'BackgroundColor', 'k',...
    % 'color',[1,1,1]);
    % box off
    % SigTimeBox(gca, xtime(1),xtime(end), [mean(confint(:,1)),  mean(confint(:,2))],'k');