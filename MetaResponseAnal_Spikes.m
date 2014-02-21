function [out_struct, hfig] = MetaResponseAnal_Spikes(input_struct)
%feed in only the highpass data you want spikes from
%don't deal with parsing baseline versus stimulus spikes in here

allfields = fieldnames(input_struct);
for ifield = 1:size(allfields,1)
    s = [allfields{ifield} ' = input_struct.' allfields{ifield} ';'];
    eval(s)
end

[spikesmat, gausstosmooth]=getspikesmat(highpassdata,spk_thresh,dt);
spiketimes = [];
for itrial=1:size(spikesmat,1)
    spiketimes{itrial}=find(spikesmat(itrial,:));
end

spkwin_size = round((0.05/dt)/2);
vm_win = [0.005, 0.002]; %this is 3 milliseconds a little before the spike;
out_struct.spk_shape = [];
out_struct.spk_vm = [];
out_struct.spk_thrsh = [];
spk_ind = 1;
for itrial=1:size(spikesmat,1)
    spks_trial = spiketimes{itrial};
    for ispike = 1:size(spks_trial,2)
        t1 = spks_trial(ispike);
        spk_wintight = [t1-20,t1+20];
        spk_win = [(t1 - spkwin_size),(t1 + spkwin_size)];
%         pad_begin = [];
%         pad_end = [];
        if spk_win(1) < 0 || spk_win(2) > size(sigdata,2)
            continue
%             pad_begin = NaN(1,0 - spk_win(1)+1);
%             spk_win(1) = 1;
%         end
%         if spk_win(2) > size(sigdata,2)
%             pad_end = NaN(1, spk_win(2) - size(sigdata,2) );
%             spk_win(2) = size(sigdata,2);
        end
        tmpshape = sigdata(itrial,spk_win(1):spk_win(2)); 
%         out_struct.spk_shape(spk_ind,:) = [pad_begin, ...
%             sigdata(itrial,spk_win(1):spk_win(2)), pad_end];
        %get peak in this window and re-center based on real spike peak
        spk_peak(spk_ind) = min(find(tmpshape == ...
            max(sigdata(itrial,spk_wintight(1):spk_wintight(2))))+spk_win(1));
        spk_win = [(spk_peak(spk_ind) - spkwin_size),(spk_peak(spk_ind) + spkwin_size)];
        spk_wintight = [spk_peak(spk_ind)-15,spk_peak(spk_ind)+15] - spk_win(1);
       
        if spk_win(1) < 0 || spk_win(2) > size(sigdata,2)
           continue
        end
        out_struct.spk_win(spk_ind,:) = spk_win;
        out_struct.spk_shape(spk_ind,:) = ...
            sigdata(itrial,spk_win(1):spk_win(2));
%         shift_x = find(out_struct.spk_shape(spk_ind,:)==max(out_struct.spk_shape(spk_ind,:)));
        
        dV_win = round(276 - (vm_win/dt)) + spk_win(1);
        out_struct.spk_vm(spk_ind) = mean(mean(sigdata(itrial, dV_win(1):dV_win(2)),1));
        dV = sigdata(itrial, dV_win(2)) - sigdata(itrial, dV_win(1));
        out_struct.spk_dv(spk_ind) = dV / (diff(dV_win)*dt);
        
        %get threshold automatically using peak and dv/dt
        dvdt = diff(out_struct.spk_shape(spk_ind,:));
        dvdt_max = max(dvdt(1,spk_wintight(1):spk_wintight(2)));
        peakind = find(dvdt == dvdt_max) ;
%         spk_peak(spk_ind) - spk_win(1);
%         dvdt_thr = 0.033 * dvdt_max;
        dvdt_thr = 0.1 * dvdt_max;
        for isamp = 1:20
            if dvdt(1,peakind - isamp) < dvdt_thr
               out_struct.spk_thrsh(spk_ind) = out_struct.spk_shape(spk_ind, (peakind - isamp + 1));
                break
            end
        end 
        spk_ind = spk_ind +1;
        
    end
end


[colvec colsize rowvec rowsize] = subplotinds(1,3);
hfig = [];
% plot spikes
if ~isempty(out_struct.spk_shape)
%     shift_x = find(mean(out_struct.spk_shape)==max(mean(out_struct.spk_shape)));
    shift_x = unique(spk_peak' - out_struct.spk_win(:,1));
    xtime = ([1:size(out_struct.spk_shape,2)]-shift_x)*dt;
    hfig(1) = figure;
    hold on
    hs1 = subplot('Position',[colvec(1),rowvec(1),colsize,rowsize]);
    line(xtime,out_struct.spk_shape')
    axis tight
    hs2 = subplot('Position',[colvec(2),rowvec(2),colsize,rowsize]);
    line(xtime,nanmean(out_struct.spk_shape)')
    axis tight
    text(-0.02,-30,['# spikes = ' num2str(size(out_struct.spk_thrsh,2))])
    hs3 = subplot('Position',[colvec(3),rowvec(3),colsize,rowsize]);
    edges = [-80 : 2: -20];
    nthrsh = histc(out_struct.spk_thrsh,edges);
    bar(edges, nthrsh/size(out_struct.spk_thrsh,2),'k')
    set(gca,'XLim',[-70,-20],'YLim',[0,1]);
    text(-60,0.8,['mean threshold = ' num2str(mean(out_struct.spk_thrsh))])
    set(hfig(1),'Position',[440   134   665   664]);
end

