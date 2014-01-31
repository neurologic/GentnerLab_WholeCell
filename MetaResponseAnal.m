r=rigdef('mac')
%%


save([r.Dir.Expt 'Analysis/Meta_10Trial_Responses_DB/exptsForICanal.mat'],...
    'exptnames','readyexpt','sigdata','vardata','basewin','allstims','repexpts',...
    'addstimcond','vmresp','base_vmresp','base_spkresp','spkresp','respDB',...
    'respdata','cb','cr','respslope','baseslope','nresid_stim','nresid_base',...
    's_base','s_stim','base_spk_vm','stim_spk_vm','stim_spk_thresh','stim_spk_thresh')

%%
load([r.Dir.Expt 'Analysis/Meta_10Trial_Responses_DB/exptsForICanal.mat'])

%%
repexpts = {
'KP_B130_131120_p1c1a.mat'
'KP_B130_131120_p1c2a.mat'
'KP_B130_131120_p1c2b.mat'
'KP_B130_131120_p1c2c.mat'
'KP_B130_131120_p1c2d.mat'
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

% trials = 5;
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

PlotRinRsCm = 0;
plotResidVm = 0;
for iexpt=1:size(repexpts,1)
    iexpt
    thisexpt=repexpts{iexpt};
    load([r.Dir.Expt thisexpt])
    
    vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
    
    table=getClampTab(expt,{'clamp',0});
    
    
    highpassdata=HighpassGeneral(vmexpt.wc.data,[],1/expt.wc.dt);
    
    outcell=PlotAndAsk(highpassdata*1000,'spk_thresh','negative');
    cellfun(@eval,outcell);
    %if spikes are not big enough... continue and dond't use this expt
    if spk_thresh < 10
        continue
    end
    
    keepsigs=reprequire(table,trials);
    allsig=table.sigsplayed;
    stimcond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
    repsigs=allsig(keepsigs);
    basetimes = expt.analysis.params.baselinewin;
    
    
    [sigon,sigoff]=GetSigTimes(expt,expt.stimcond,1);
    
    %%%%%%% need to go through and give non-dBexpt stimuli a dummy db so
    %%%%%%% can include them in these analyses but group them appropriately
    %%%%%% can i even go back and measure what dB they are?
%     [dbstimcond,dblevels]=getDBstimcond(vmexpt);
%     if isempty(dbstimcond)
%         
%         notdbexpt(iexpt) = iexpt;
%         continue
%     end
    
    %set up input_struct with data to pass analysis functions
    input_struct.basetimes = basetimes;
    input_struct.sigon = sigon;
    input_struct.sigoff = sigoff;
    input_struct.spk_thresh = spk_thresh;

    %get input and access resistance from the current step preceding trial
    [out_struct, hfig] = MetaResponseAnal_RsRin(vmexpt,input_struct);
    allfields = fieldnames(out_struct);
    for ifield = 1:size(allfields,1)
        s = [allfields{ifield} '(iexpt) = out_struct.' allfields{ifield} ';'];
        eval(s)
    end
    if PlotRinRsCm == 1
        foldername = '/Users/kperks/GitHub/Data_Mat/Analysis/Meta_10Trial_Responses_DB/RinRsCm/';
        saveas(hfig,[foldername expt.name 'Stepdata_smartCorrect.fig'])
        saveas(hfig,[foldername expt.name 'Stepdata_smartCorrect.png'])
    end
    
    close(hfig)
    

    physiol = input('are these estimates physiological? 1:yes 0:no')
    if physiol == 0
       Rin(iexpt) = NaN;
       Rs(iexpt) = NaN;
       Cm(iexpt) = NaN;
    end
    close(hfig)
    
    %go through each relevant stimulus for this cell...
    % do this set of analysis for each type of dB signal
    for istimcond=1:size(dbstimcond,2)
        exptID(stimind) = iexpt;
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
        
        %%%%%%%%%%%%%         for now just use a single dB stim set at a time
        dbind = find(thisdb==signalDB);
        if isempty(dbind)
            continue
        end
        
        sigexpt = filtesweeps(vmexpt,0,'wavnames',...
            thiscond(dbind).wavnames);
        sigdata = sigexpt.wc.data*1000;
        sigdata_filt = medfilt1(sigdata,200,[],2);
         
        %set up structure of data to send to all analyses
        input_struct.sigdata = sigdata;
        input_struct.sigdata_filt = sigdata_filt;
         
        %get power spectrum
        out_struct = MetaResponseAnal_PDS(expt, input_struct)      
        allfields = fieldnames(out_struct);
        for ifield = 1:size(allfields,1)
            s = [allfields{ifield} ' = out_struct.' allfields{ifield} ';'];
            eval(s)
        end
        
        %get residuals of Vm distributions for baseline and stimulus period
        MetaResponseAnal_VmResiduals(expt,input_struct)
        allfields = fieldnames(out_struct);
        for ifield = 1:size(allfields,1)
            s = [allfields{ifield} ' = out_struct.' allfields{ifield} ';'];
            eval(s)
        end
        
        % get spike shapes and spike threshold and actual spike peak times to get vm during spike
        MetaResponseAnal_Spikes(expt,input_struct)
        allfields = fieldnames(out_struct);
        for ifield = 1:size(allfields,1)
            s = [allfields{ifield} ' = out_struct.' allfields{ifield} ';'];
            eval(s)
        end
        
        % get vm response windows
        MetaResponseAnal_VmResponseWin(expt,input_struct)
        allfields = fieldnames(out_struct);
        for ifield = 1:size(allfields,1)
            s = [allfields{ifield} ' = out_struct.' allfields{ifield} ';'];
            eval(s)
        end
        
        stimind = stimind + 1;
    end
    
end
%%
%plot population distributions for RinRsCm
RinBins = [0:25:400];
RinCT = histc(Rin,RinBins)/size(Rin,2);
figure;
line(RinBins,RinCT);
ylabel('Rin')
[x,p] = empcdf (Rin);
figure;
stairs(x,p)

RsBins = [0:10:100];
RsCT = histc(Rs,RsBins)/size(Rs,2);
figure;
line(RsBins,RsCT);
ylabel('Rs')
[x,p] = empcdf (Rs);
figure;
stairs(x,p)

CmBins = [0:20:800];
CmCT = histc(Cm,CmBins)/size(Cm,2);
figure;
line(CmBins,CmCT);
ylabel('Cm')
[x,p] = empcdf (Cm);
figure;
stairs(x,p)

figure;scatter(response_vm,response_spk,50,'k','fill')
ylabel('average spike rate per response')
xlabel('average membrane potential per "up_win" response')
title ('16 cells; 210 "up" responses; 40dB stimuli only')

%mean distribution of Vm residuals
figure;
hold on
stairs(resid_edges,mean(nresid_base),'color','b','LineWidth',3)
stairs(resid_edges,mean(nresid_stim),'color','r','LineWidth',3)
axis tight
ylims = get(gca,'YLim');
text(resid_edges(10),ylims(2)-0.05,{['mean base skew = ' num2str(mean(skew_base))];...
    ['mean stim skew = ' num2str(mean(skew_stim))]},'HorizontalAlignment','center');
title([num2str(stimind-1) ' cell:signal pairs    stimulus at ' ...
    num2str(signalDB) 'dBSPL'])

% spikethresholds for base and stim
valind = find(stim_spk_thresh~=0);
stim_spk_thresh = stim_spk_thresh(valind);

valind = find(base_spk_thresh~=0);
base_spk_thresh = base_spk_thresh(valind);

figure;
hold on
scatter(repmat(1,size(base_spk_thresh,2),1),base_spk_thresh,50,'b','fill');
scatter(repmat(2,size(stim_spk_thresh,2),1),stim_spk_thresh,50,'r','fill');

set(gca,'XTick',[1,2],'XLim',[0,3],'XTickLabel',{'baseline','stimulus'})
ylabel('spike  (mV)')

%using residuals is a really easy way to "highpass"

% %p(spk|Vm)
% figure;
% hold on
% scatter(x_bins,nanmean(P_spk_base,1),100,'b','fill')
% scatter(x_bins,nanmean(P_spk_resp,1),100,'r','fill')
% xlabel('Vm')
% ylabel('normalized p(spike|Vm)')
% set(gca,'XTick', x_bins, 'XTickLabel', x_bins)
% title([num2str(stimind-1) ' cell:signal pairs    stimulus at ' ...
%     num2str(signalDB) 'dBSPL'])
