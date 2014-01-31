function MetaResponseAnal_Spikes(expt,input_struct)

allfields = fieldnames(input_struct);
for ifield = 1:size(allfields,1)
   s = [allfields{ifield} ' = input_struct.' allfields{ifield} ';'];
   eval(s)
end

highpassdata=HighpassGeneral(sigdata,[],1/expt.wc.dt);
[spikesmat, gausstosmooth]=getspikesmat(highpassdata,spk_thresh,expt);
spiketimes = [];
for itrial=1:size(spikesmat,1)
    spiketimes{itrial}=find(spikesmat(itrial,:));
end

spkwin_size = round((0.05/expt.wc.dt)/2);
basespk_ind = 1;
stimspk_ind = 1;
base_spk_vm{stimind} = [];
stim_spk_vm{stimind} = [];

for itrial=1:size(spikesmat,1)
    spks_trial = spiketimes{itrial};
    for ispike = 1:size(spks_trial,2)
        t1 = spks_trial(ispike);
        %if the spike is during baseline
        if t1 > basetimes(1) && t1 < sigon
            spk_win = [(t1 - spkwin_size),(t1 + spkwin_size)];
            base_spk_vm{stimind}(basespk_ind,:) = sigdata(itrial,spk_win(1):spk_win(2));
            %get peak in this window and re-center based on real spike
            %peak
            spk_peak = find(base_spk_vm{stimind}(basespk_ind,:) == max(base_spk_vm{stimind}(basespk_ind,:)))+spk_win(1);
            spk_win = [(spk_peak - spkwin_size),(spk_peak + spkwin_size)];
            base_spk_vm{stimind}(basespk_ind,:) = sigdata(itrial,spk_win(1):spk_win(2));
            basespk_ind = basespk_ind +1;
        end
        %if the spike is during the stimulus (try + 200msec?
        % given what i know about the long-lasting recovery and offset resp?...
        if t1 > sigon && t1 < sigoff
            spk_win = [(t1 - spkwin_size),(t1 + spkwin_size)];
            stim_spk_vm{stimind}(stimspk_ind,:) = sigdata(itrial,spk_win(1):spk_win(2));
            %get peak in this window and re-center based on real spike
            %peak
            spk_peak = min(find(stim_spk_vm{stimind}(stimspk_ind,:) == max(stim_spk_vm{stimind}(stimspk_ind,:)))+spk_win(1));
            spk_win = [(spk_peak - spkwin_size),(spk_peak + spkwin_size)];
            stim_spk_vm{stimind}(stimspk_ind,:) = sigdata(itrial,spk_win(1):spk_win(2));
            stimspk_ind = stimspk_ind +1;
        end
    end
end

if plotSpikeShape == 1;
    [colvec colsize rowvec rowsize] = subplotinds(1,2);
    hfig = [];
    % plot baseline spikes
    if ~isempty(base_spk_vm{stimind})
        shift_x = find(mean(base_spk_vm{stimind})==max(mean(base_spk_vm{stimind})));
        xtime = ([1:size(base_spk_vm{stimind},2)]-shift_x)*expt.wc.dt;
        hfig(1) = figure;
        hold on
        hs1 = subplot('Position',[colvec(1),rowvec(1),colsize,rowsize]);
        line(xtime,base_spk_vm{stimind}')
        axis tight
        title(['expt.name : stimulus #' num2str(istimcond)...
            '    # baseline spikes = ' num2str(size(base_spk_vm{stimind},1))])
        hs2 = subplot('Position',[colvec(2),rowvec(2),colsize,rowsize]);
        line(xtime,mean(base_spk_vm{stimind})')
        axis tight
        text(-.02,0,['expt.name : stimulus #' num2str(istimcond)...
            '    # baseline spkies = ' num2str(size(base_spk_vm{stimind},1))])
        set(hfig(1),'Position',[440   134   665   664]);
        saveas(hfig(1),[r.Dir.Expt 'Analysis/Meta_10Trial_Responses_DB/SpikeShapes/' ...
            expt.name '_allBaseSpk.png'])
        saveas(hfig(1),[r.Dir.Expt 'Analysis/Meta_10Trial_Responses_DB/SpikeShapes/' ...
            expt.name '_allBaseSpk.fig'])
        base_spk_mean(stimind,:) = mean(base_spk_vm{stimind});
        spk_thresh(stimind,1) = input('spike threshold?');
    end
    
    %plot stimulus spikes
    if ~isempty(stim_spk_vm{stimind})
        shift_x = find(mean(stim_spk_vm{stimind})==max(mean(stim_spk_vm{stimind})));
        xtime = ([1:size(stim_spk_vm{stimind},2)]-shift_x)*expt.wc.dt;
        hfig(2) = figure;
        hold on
        hs1 = subplot('Position',[colvec(1),rowvec(1),colsize,rowsize]);
        line(xtime,stim_spk_vm{stimind}')
        axis tight
        title(['expt.name : stimulus #' num2str(istimcond)...
            '    # stim spkies = ' num2str(size(stim_spk_vm{stimind},1))])
        hs2 = subplot('Position',[colvec(2),rowvec(2),colsize,rowsize]);
        line(xtime,mean(stim_spk_vm{stimind})')
        axis tight
        set(hfig(2),'Position',[440   134   665   664]);
        saveas(hfig(2),[r.Dir.Expt 'Analysis/Meta_10Trial_Responses_DB/SpikeShapes/' ...
            expt.name '_allStimSpk.png'])
        saveas(hfig(2),[r.Dir.Expt 'Analysis/Meta_10Trial_Responses_DB/SpikeShapes/' ...
            expt.name '_allStimSpk.fig'])
        stim_spk_mean(stimind,:) = mean(stim_spk_vm{stimind});
        spk_thresh(stimind,2) = input('spike thresh ?');
    end
end
spk_bad(stimind) = input('spikes bad (0) or good (1)?');
for ifig = 1:size(hfig,2)
    if hfig(ifig)~=0
        close(hfig(ifig))
    end
end
hfig = [];