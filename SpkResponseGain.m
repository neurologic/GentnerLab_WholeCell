% get psth from batchbatchPSTHionto.m
%
% basedata = base_data_All;
% gzdata = gz_data_All;

%in JVT_getPSTH, I have it so that it calculates the psth only for the
%stimulus period... 
function Result = SpkResponseGain(basedata,gzdata,sigon,bps)
Result.baseRate_resp = [];
Result.gzRate_resp = [];
Result.Gain_resp = [];
Result.baseRate_pre = [];
Result.gzRate_pre = [];
Result.Gain_pre = [];

plotall = 1;
respind = 1;
for istim = 1:size(basedata,2)
    thisbase = basedata(istim).spikes;
    thisgz = gzdata(istim).spikes;
    usebase = thisgz;
    %return the Result struct with:
    %     gain of each response estimated just by dividing
    %     estimated gain by fit line when more than one response per cell
    
    %get response windows off of basedata
    % highpassdata = basedata;
    % spk_thresh = 0.5;
    dt = 1/bps;
    binwidth = 0.05; %seconds
    % [out_struct, hfig] = MetaResponseAnal_Spikes(input_struct);
    %oops. wanted to use that function to keep things standardized, but that
    %funtion is for working with highpass and lowpass filtered versions of
    %intracellular recorded data for getting stuff like spike shape and
    %threshold and p(spk|vm)
    % [spikesmat, gausstosmooth]=getspikesmat(highpassdata,spk_thresh,dt);
    %don't need to get spikesmat because basedata from psth script is already
    %spkvec!
    % spkvec = zeros(size(spikesmat,1),size(spikesmat,2));
    % spkvec(find(spikesmat)) = 1;
    
    for itrial = 1:size(usebase,1)
        trialisi(itrial) = median(diff(find(usebase(itrial,:) ==1)));
    end
    sampisi = round(mean(trialisi));
    if isnan(sampisi) || sampisi == 0 || sampisi >= sigon 
        sampisi = 0;
        gausstosmooth = fspecial('gaussian', [200*8,1],200);
        gausstosmooth = gausstosmooth/max(gausstosmooth);
    end
    if ~isnan(sampisi) && sampisi > 0;
        %     gausstosmooth = fspecial('gaussian', [allisi*4,1],allisi/2);
        gausstosmooth = fspecial('gaussian', [sampisi*8,1],sampisi);
        gausstosmooth = gausstosmooth/max(gausstosmooth);
    end
    
    smoothspk = [];
    for itrial = 1:size(usebase,1)
        smoothspk(itrial,:) = conv(usebase(itrial,:),gausstosmooth,'same');
    end
    %                 input_struct.sigdata_filt = smoothspk;
    
    conf_p = 95;
    windowsize = binwidth / dt; %samples
    binsize = 10; %smaples
    bin_p = 0.85;
    
%     quantile(mean(smoothspk(:,1:sigon-sampisi)),conf_p);
    confint_spk = getCDFconf (mean(smoothspk(:,1:sigon-sampisi)),conf_p);
    
    [up_inds_spk, low_inds_spk] = WindowResponse(mean(smoothspk(:,sigon:end)),...
        confint_spk, windowsize, binsize, bin_p);
    up_win_spk = getWindowEdges (up_inds_spk, 1, 1)+sigon;
    low_win_spk = getWindowEdges (low_inds_spk, 1, 1)+sigon;
    
    if plotall ==1
        spkrespfig = figure;
        hold on
        subplot(2,1,1)
        hold on
        for itrial = 1:size(smoothspk,1)
            thesespks = find(thisbase(itrial,:))*dt;
            for ispike = 1:size(thesespks,2)
                line([thesespks(ispike),thesespks(ispike)],...
                    [1*itrial,(1*itrial)+1],'color','k')
            end
        end
        set(gca, 'YLim',[1,size(smoothspk,1)+1],'XLim',[1*dt,size(smoothspk,2)*dt])
        
        SigTimeBox(gca, (sigon)*dt,size(thisbase,2)*dt, get(gca,'YLim'),[0.5 0.5 0.5]);
        for iresp=1:size(up_win_spk,2)
            SigTimeBox(gca, up_win_spk(1,iresp)*dt, ...
                up_win_spk(2,iresp)*dt, get(gca,'YLim'),'r');
        end
        for inhib=1:size(low_win_spk,2)
            SigTimeBox(gca, low_win_spk(1,inhib)*dt, ...
                low_win_spk(2,inhib)*dt, get(gca,'YLim'),'b');
        end
        axis tight
        set(gca,'TickDir','out')
        box off
        
        subplot(2,1,2)
        hold on
        for itrial = 1:size(smoothspk,1)
            thesespks = find(thisgz(itrial,:))*dt;
            for ispike = 1:size(thesespks,2)
                line([thesespks(ispike),thesespks(ispike)],...
                    [1*itrial,(1*itrial)+1],'color','k')
            end
        end
        set(gca, 'YLim',[1,size(smoothspk,1)+1],'XLim',[1*dt,size(smoothspk,2)*dt])
        
        SigTimeBox(gca, (sigon)*dt,size(thisbase,2)*dt, get(gca,'YLim'),[0.5 0.5 0.5]);
        for iresp=1:size(up_win_spk,2)
            SigTimeBox(gca, up_win_spk(1,iresp)*dt, ...
                up_win_spk(2,iresp)*dt, get(gca,'YLim'),'r');
        end
        for inhib=1:size(low_win_spk,2)
            SigTimeBox(gca, low_win_spk(1,inhib)*dt, ...
                low_win_spk(2,inhib)*dt, get(gca,'YLim'),'b');
        end
        axis tight
        set(gca,'TickDir','out')
        box off
        
        set(spkrespfig,'Position',[212         523        1168         283])
    end
    
    
    numbasetrials = size(thisbase,1);
    numgztrials = size(thisgz,1);
    for iresp=1:size(up_win_spk,2)
        Tstart = up_win_spk(1,iresp);
        Tstop = up_win_spk(2,iresp);
        Result.baseRate_resp(respind) = size(find(thisbase(:,Tstart:Tstop)),1)/numbasetrials/((Tstop-Tstart)*dt);
        Result.gzRate_resp(respind) = size(find(thisgz(:,Tstart:Tstop)),1)/numbasetrials/((Tstop-Tstart)*dt);
        Result.Gain_resp(respind) = Result.gzRate_resp(respind) / Result.baseRate_resp(respind);
        % spike time distribution... for each response, all times are
        % normalized to max time so that histograms are in quantiles of the
        % response
        basespkT = [];
        gzspkT = [];
        edges = [0:0.05:1];
        for itrial = 1:size(thisbase,1)
            tmpbasespksT = find(thisbase(itrial,Tstart:Tstop))*dt;
            tmpbasespksT  = tmpbasespksT * (1/ ((Tstop-Tstart)*dt));
            basespkT = [basespkT, tmpbasespksT];
            
            tmpgzspksT = find(thisgz(itrial,Tstart:Tstop))*dt;
            tmpgzspksT = tmpgzspksT * (1/ ((Tstop-Tstart)*dt));
            gzspkT = [gzspkT, tmpgzspksT];
        end
%         Result.basespkHist {istim}(iresp,:) = histc(basespksT,edges)/size(basespksT,2);
        Result.basespkT {respind} = basespkT;
%         Result.gzspkHist {istim}(iresp,:) = histc(gzspksT,edges)/size(gzspksT,2);
        Result.gzspkT {respind} = gzspkT;
        respind = respind +1;
    end
    Result.baseRate_pre(istim) = size(find(thisbase(:,1:sigon)),1)/numbasetrials/(sigon*dt);
    Result.gzRate_pre(istim) = size(find(thisgz(:,1:sigon)),1)/numbasetrials/(sigon*dt);
    Result.Gain_pre(istim) = Result.gzRate_pre(istim) / Result.baseRate_pre(istim);
    stimnum = istim
end
%don't do inhib or no spikes "yet" (if ever)
%                 for inhib=1:size(low_win_spk,2)
%                     SigTimeBox(gca, low_win_spk(1,inhib)*expt.wc.dt, ...
%                         low_win_spk(2,inhib)*expt.wc.dt, get(gca,'YLim'),'b');
%                 end

%
%plot without patches for eps
%                 spkrespfig = figure;
%                 hold on
%                 for itrial = 1:size(smoothspk,1)
%                     thesespks = find(spkvec(itrial,:))*expt.wc.dt
%                    for ispike = 1:size(thesespks,2)
%                        line([thesespks(ispike),thesespks(ispike)],[1*itrial,(1*itrial)+1])
%                    end
%                 end
%                 set(gca, 'YLim',[-0.5,size(smoothspk,1)+1],'XLim',[1*expt.wc.dt,size(smoothspk,2)*expt.wc.dt])
%
%                 line([(sigon)*expt.wc.dt,sigoff*expt.wc.dt], [1,1],'color','k');
%                 for iresp=1:size(up_win_spk,2)
%                     line([up_win_spk(1,iresp)*expt.wc.dt, ...
%                         up_win_spk(2,iresp)*expt.wc.dt], [0.5,0.5],'color','r');
%                 end
%                 for inhib=1:size(low_win_spk,2)
%                     line([low_win_spk(1,inhib)*expt.wc.dt, ...
%                         low_win_spk(2,inhib)*expt.wc.dt], [0,0],'color','b');
%                 end
%                   box off
%                   set(gca,'TickDir','out')
%                 set(spkrespfig,'Position',[212         523        1168         283])
%                 title([expt.name stimcond(istim).wavnames],'Interpreter','none')
% %
end

