

fs = 44100;
max_lag_c = [];
max_lag_t = [];
min_lag_c = [];
min_lag_t = [];
allresp = [];
allrms = [];
verboseresp = [];

allind = 1;
hfig = figure
hold on
for iunq = 1:size(unq_expts,2)
    this_cell = unq_expts{iunq}
    sigind = 1;
    %do an averate across expts and stims for each unique cell
    for iexpt=1:size(repexpts,1)
        if ~isempty(regexp(repexpts{iexpt},unq_expts{iunq}))
            thisexpt=repexpts{iexpt};
            load([r.Dir.Expt thisexpt])
            vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
            table=getClampTab(expt,{'clamp',0});
            keepsigs=reprequire(table,trials);
            thiscond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
    
            for isig = 1:size(thiscond,2)
                %skip warped shortened tempo stims
                if ~isempty(regexp(thiscond(isig).wavnames,'ws'))
                    message = 'is warped stim'
                    continue
                end
                %%%%%%%%%%%%i guess i high-pass fft'ed these at some point?
                %         A = {'one','two','twenty-two','One','two'};
                sigexpt = filtesweeps(vmexpt,0,'wavnames',thiscond(isig).wavnames);
                verboseresp = medfilt1(sigexpt.wc.data,200,[],2);
                [sigon,sigoff]=GetSigTimes(sigexpt,thiscond,isig);
                baselinewin = expt.analysis.params.baselinewin;
                
                prevar = var(verboseresp(:,baselinewin(1):sigon));
                stimvar = var(verboseresp(:,sigon:sigoff));
                %if the variance distribution during stim is different than
                %pre-stim... then that counts as a response
                if kstest2(prevar,stimvar)==1
                    y = thiscond(isig).wavs;
                    %             y = wavread(fullfile(r.Dir.Stims,thiscond(isig).wavnames));
                    %             filtwav = fftFilter(y,fs,50,2);
                    %             y = filtwav;
                    dcoff = (mean(y,1));
                    nodc = y-dcoff;
                    tmprms = sqrt((nodc.^2));
                    olddt=1/fs;
                    bin=round(1/olddt/round(1/expt.wc.dt));
                    newdt=olddt*bin;
                    tmprms=tmprms(1:bin:end,:);
                     tmprms = tmprms(1:minlen) - min(tmprms(1:minlen));
                     
                     tmprms = tmprms / max(tmprms);
                      sampsmooth = round(0.002/newdt);
                   gausstosmooth = fspecial('gaussian', [sampsmooth*8,1],sampsmooth/4);
                   gausstosmooth = gausstosmooth / max(gausstosmooth);
                   tmprms = conv(tmprms,gausstosmooth,'same');
                   tmprms = tmprms / max(tmprms);
                    allrms(allind,:) = tmprms';
                    
                    tmpresp = mean(verboseresp(:,sigon:sigon+minlen));
                    tmpresp = tmpresp-min(tmpresp);
                    tmpresp = tmpresp / max(tmpresp);
                    allresp(allind,:) = tmpresp;
                    
                    tmpvar = var(verboseresp(:,sigon:sigon+minlen));
                    
                    [XC,LAGS] = xcorr(tmpresp,tmprms');
                    xc = XC(1,find(LAGS > 200));
                    lags = LAGS(1,find(LAGS > 200))*expt.wc.dt;
                    max_lag_c{iunq}(sigind) = xc(min(find(xc == max(xc))));
                    max_lag_t{iunq}(sigind) = lags(min(find(xc == max(xc))));
                    min_lag_c{iunq}(sigind) = xc(min(find(xc == min(xc))));
                    min_lag_t{iunq}(sigind) = lags(min(find(xc == min(xc))));
                    sigind = sigind + 1;
                    allind = allind + 1;
                end
            end
            
        end
    end
    
    subplot(ceil(size(unq_expts,2)/4),4,iunq)
    hold on
    scatter(max_lag_t{iunq},max_lag_c{iunq},50,'r','fill')
%     scatter(min_lag_t{iunq},min_lag_c{iunq},50,'b','fill')
%     
    scatter(mean(max_lag_t{iunq}),mean(max_lag_c{iunq}),100,'k','fill')
    qc = quantile(max_lag_c{iunq},[0.15,0.85]);
    qt = quantile(max_lag_t{iunq},[0.15,0.85]);
    line([mean(max_lag_t{iunq}),mean(max_lag_t{iunq})],qc,...
        'color','k','LineWidth',2)
    line(qt,[mean(max_lag_c{iunq}),mean(max_lag_c{iunq})],...
        'color','k','LineWidth',2)
    set(gca,'XLim',[0,0.6],'YLim',[0,2250])
    title(this_cell,'Interpreter','none','FontSize',14)
    %     legend('max corr','min corr')
    
end

for iunq = 1:size(max_lag_c,2)
    isnorm_maxc(iunq) = kstest(max_lag_c{iunq});
    isnorm_maxl(iunq) = kstest(max_lag_t{iunq});
end
%enough are not norm that will plot bars as 15% and 85% quantiles

figure;
hold on
for iunq = 1:size(max_lag_c,2)
    scatter(mean(max_lag_t{iunq}),mean(max_lag_c{iunq}),100,'k','fill')
    qc = quantile(max_lag_c{iunq},[0.15,0.85]);
    qt = quantile(max_lag_t{iunq},[0.15,0.85]);
    line([mean(max_lag_t{iunq}),mean(max_lag_t{iunq})],qc,...
        'color','k','LineWidth',2)
    line(qt,[mean(max_lag_c{iunq}),mean(max_lag_c{iunq})],...
        'color','k','LineWidth',2)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now randomize across responses and stims to get "noise" bounds on estimate
%of correlation and lags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randmax_lag_c = [];
randmax_lag_t = [];
randmin_lag_c = [];
randmin_lag_t = [];
randmax_c_Mn = [];
randmax_t_Mn = [];
monind = [1:size(allresp,1)];

numrandomize = 5;
for iunq = 1:size(unq_expts,2)
    numsigs = size(max_lag_c{iunq},2);
%     randind = randsample(size(allresp,1),size(allresp,1))';
    for irand = 1:numrandomize
    originvec = [1:size(allresp,1)];
    randindA = randsample(size(originvec,2),numsigs);
    originvec(randindA) = [];
    originvec = find(originvec);
    randindB = randsample(size(originvec,2),numsigs);
    
    for isig = 1:numsigs
        thisrms = allrms(randindA(isig),:);
        thisresp = allresp(randindB(isig),:);
        
        [XC,LAGS] = xcorr(thisresp,thisrms);
        xc = XC(1,find(LAGS > 200));
        lags = LAGS(1,find(LAGS > 200))*expt.wc.dt;
        randmax_lag_c(isig) = xc(min(find(xc == max(xc))));
        randmax_lag_t(isig) = lags(min(find(xc == max(xc))));
        randmin_lag_c(isig) = xc(min(find(xc == min(xc))));
        randmin_lag_t(isig) = lags(min(find(xc == min(xc))));
        
    end
    randmax_c_Mn{iunq}(irand) = mean(randmax_lag_c);
    randmax_t_Mn{iunq}(irand) = mean(randmax_lag_t);
    end
end

hfig = figure
hold on
for iunq = 1:size(unq_expts,2)
    this_cell = unq_expts{iunq}  
    subplot(ceil(size(unq_expts,2)/4),4,iunq)
    hold on
    scatter(max_lag_t{iunq},max_lag_c{iunq},50,'r','fill')
%     scatter(min_lag_t{iunq},min_lag_c{iunq},50,'b','fill')
%     
    scatter(mean(max_lag_t{iunq}),mean(max_lag_c{iunq}),100,'k','fill')
    qc = quantile(max_lag_c{iunq},[0.15,0.85]);
    qt = quantile(max_lag_t{iunq},[0.15,0.85]);
    line([mean(max_lag_t{iunq}),mean(max_lag_t{iunq})],qc,...
        'color','k','LineWidth',2)
    line(qt,[mean(max_lag_c{iunq}),mean(max_lag_c{iunq})],...
        'color','k','LineWidth',2)
    set(gca,'XLim',[0,0.6],'YLim',[0,2250])
    title(this_cell,'Interpreter','none','FontSize',14)
    
    scatter(randmax_t_Mn{iunq},randmax_c_Mn{iunq},50,'b','fill')
end
%%%%%%%%%%%
% figure;hold on
scatter(mean(randmax_lag_t),mean(randmax_lag_c),100,'r','fill')
    qc = quantile(randmax_lag_c,[0.15,0.85]);
    qt = quantile(randmax_lag_t,[0.15,0.85]);
    line([mean(randmax_lag_t),mean(randmax_lag_t)],qc,...
        'color','r','LineWidth',2)
    line(qt,[mean(randmax_lag_c),mean(randmax_lag_c)],...
        'color','r','LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%
%didn't look "good" at first...
%but don't just pay attention to the mean!!!
%notice the symmetry of the CI in the random condition
%versus the assymetry for the cells on the same scale
%%%%%%%%%%%%%%%%%%%%%%

%% dump the following way to do it because don't know reasonable way to address the correlation in this way...
%trying the analysis that is now above this one to look at lags and
%correlation
sig_wavname = [];
sig_wav = [];
rmswav = [];
tmpnames = [];
sigind = 1;
fs = 44100;
response_offset = 0;%0.03; %seconds

unqRMS = [];
unqRESP = [];
unqVAR = [];

u_mn_erms = [];
u_mn_eresp = [];
u_mn_evar = [];

for iunq = 1:size(unq_expts,2)
    this_cell = unq_expts{iunq}
    
    sortRMS = [];
    sortRESP = [];
    sortVAR = [];
    mn_erms = [];
    mn_eresp = [];
    mn_evar = [];
    exptind = 1;
    for iexpt=1:size(repexpts,1)
        if ~isempty(regexp(repexpts{iexpt},unq_expts{iunq}))
            thisexpt=repexpts{iexpt};
            load([r.Dir.Expt thisexpt])
            vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
            table=getClampTab(expt,{'clamp',0});
            keepsigs=reprequire(table,trials);
            thiscond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
            
            %do an average for each expt so that they stay separate
            
            exptrms = [];
            exptresp = [];
            exptvar = [];
            
            sigind = 1;
            for isig = 1:size(thiscond,2)
                %skip warped shortened tempo stims
                if ~isempty(regexp(thiscond(isig).wavnames,'ws'))
                    message = 'is warped stim'
                    continue
                end
                %%%%%%%%%%%%i guess i high-pass fft'ed these at some point?
                %         A = {'one','two','twenty-two','One','two'};
                sigexpt = filtesweeps(vmexpt,0,'wavnames',thiscond(isig).wavnames);
                allresp = medfilt1(sigexpt.wc.data,200,[],2);
                [sigon,sigoff]=GetSigTimes(sigexpt,thiscond,isig);
                baselinewin = expt.analysis.params.baselinewin;
                
                prevar = var(allresp(:,baselinewin(1):sigon));
                stimvar = var(allresp(:,sigon:sigoff));
                %if the variance distribution during stim is different than
                %pre-stim... then that counts as a response
                if kstest2(prevar,stimvar)==1
                    y = thiscond(isig).wavs;
                    %             y = wavread(fullfile(r.Dir.Stims,thiscond(isig).wavnames));
                    %             filtwav = fftFilter(y,fs,50,2);
                    %             y = filtwav;
                    dcoff = (mean(y));
                    nodc = y-dcoff;
                    tmprms = sqrt((nodc.^2));
                    tmprms = tmprms - min(tmprms);
                    olddt=1/fs;
                    bin=round(1/olddt/round(1/expt.wc.dt));
                    newdt=olddt*bin;
                    tmprms=tmprms(1:bin:end,:);
                    tmprms = medfilt1(tmprms(1:minlen,1),200,[],1);
                    exptrms(sigind,:) = tmprms';
                    % tmprms = tmprms ./ max(tmprms);
                    
                    offsize = round(response_offset/expt.wc.dt);
                    %             bins = [1:offsize:(size(tmprms,1))+offsize];
                    
                    tmpresp = mean(allresp(:,sigon+offsize:sigon+minlen+offsize));
                    tmpresp = tmpresp - min(tmpresp);
                    exptresp(sigind,:) = tmpresp;
                    %             tmpresp = tmpresp ./ max(tmpresp);
                    
                    tmpvar = var(allresp(:,sigon+offsize:sigon+minlen+offsize));
                    tmpvar = tmpvar - min(tmpvar);
                    exptvar(sigind,:) = tmpvar;
                    %             tmpvar = tmpvar ./ max(tmpvar);
                    
                    sigind = sigind + 1;
                end
            end
            
            mn_erms(exptind,:) = mean(exptrms) ./ max(mean(exptrms));
            mn_eresp(exptind,:) = mean(exptresp) ./ max(mean(exptresp));
            mn_evar(exptind,:) = mean(exptvar) ./ max(mean(exptvar));
            
            edges = [0:0.1:1];
            figure;hold on
            subplot(2,1,1)
            hold on
            title(expt.name,'Interpreter','none')
            xtime = [1:size(mn_erms,2)]*expt.wc.dt;
            plot([1:size(mn_erms,2)]*expt.wc.dt,mn_erms(exptind,:),'color','k')
            plot([1:size(mn_eresp,2)]*expt.wc.dt,mn_eresp(exptind,:),'color','r')
            plot([1:size(mn_evar,2)]*expt.wc.dt,mn_evar(exptind,:),'color','b')
            axis tight
            
            [nrms,bini] = histc(mn_erms(exptind,:),edges);
            [nresp] = histc(mn_eresp(exptind,:),edges);
            [nvar] = histc(mn_evar(exptind,:),edges);
            subplot(2,1,2)
            hold on
            plot(edges,nrms/size(mn_erms,2),'color','k')
            plot(edges,nresp/size(mn_eresp,2),'color','r')
            plot(edges,nvar/size(mn_evar,2),'color','b')
            legend('rms','resp','var')
            
          
          exptind = exptind + 1;
        end
    end
    
    u_mn_erms(iunq,:) = normdata(mean(mn_erms,1));
    u_mn_eresp(iunq,:) = normdata(mean(mn_eresp,1));
    u_mn_evar(iunq,:) = normdata(mean(mn_evar,1));
    
    [nrms,bini] = histc(u_mn_erms(iunq,:),edges);
    for ibin = 1:size(edges,2)
        bininds = find(bini == ibin);
        unqRMS (ibin,iunq) = mean(u_mn_erms(iunq,bininds));
        [nresp] = histc(u_mn_eresp(iunq,bininds),edges);
        [nvar] = histc(u_mn_evar(iunq,bininds),edges);
        unqRESP (ibin,:,iunq) = nresp./max(size(u_mn_eresp(iunq,bininds)));
        unqVAR (ibin,:,iunq) = nvar./max(size(u_mn_evar(iunq,bininds)));
    end
  
    %     hfig(exptind) = figure;
    %     imagesc(mean(unqRMS,2),edges,mean(unqRESP,3));
    
end

%resp anti-correlates with rms??
hfig = figure;
subplot(2,1,1)
imagesc(mean(unqRMS,2),flipud(edges),mean(unqRESP,3));
subplot(2,1,2)
scatter(mean(unqRMS,2),mean(median(unqRESP,3),2),50,'k','fill')
title('resp anti-correlates with rms??','FontSize',14)

%var anti-correlates with rms??
hfig = figure;
subplot(2,1,1)
imagesc(mean(unqRMS,2),flipud(edges),mean(unqVAR,3));
subplot(2,1,2)
scatter(mean(unqRMS,2),mean(median(unqVAR,3),2),50,'k','fill')
title('var anti-correlates with rms??','FontSize',14)

%var correlates with resp
hfig = figure;
subplot(2,1,1)
imagesc(mean(mean(unqRESP,3),1),flipud(edges),mean(unqVAR,3));
subplot(2,1,2)
scatter(mean(mean(unqRESP,3),1),mean(median(unqVAR,3),2),50,'k','fill')
title('var correlates with resp?','FontSize',14)


figure;hold on
subplot(2,1,1)
hold on
title('all unique experiemnts averaged','FontSize',14,'Interpreter','none')
xtime = [1:size(u_mn_erms,2)]*expt.wc.dt;
plot([1:size(u_mn_erms,2)]*expt.wc.dt,normdata(mean(u_mn_erms)),'color','k')
plot([1:size(u_mn_eresp,2)]*expt.wc.dt,normdata(mean(u_mn_eresp)),'color','r')
plot([1:size(u_mn_evar,2)]*expt.wc.dt,normdata(mean(u_mn_evar)),'color','b')
axis tight
edges = [0:0.1:1];
[nrms,bini] = histc(normdata(mean(u_mn_erms)),edges);
[nresp] = histc(normdata(mean(u_mn_eresp)),edges);
[nvar] = histc(normdata(mean(u_mn_evar)),edges);
subplot(2,1,2)
hold on
plot(edges,nrms/size(u_mn_erms,2),'color','k')
plot(edges,nresp/size(u_mn_eresp,2),'color','r')
plot(edges,nvar/size(u_mn_evar,2),'color','b')
legend('rms','resp','var')


%
% edges = [0:0.1:1];
% [nrms,bini] = histc(meta_rms,edges);
% sortRESP = [];
% sortVAR = [];
%
% for ibin = 1:size(edges,2)
%     bininds = find(bini == ibin);
%     sortRMS (ibin,sigind) = mean(meta_rms(bininds));
%     [nresp] = histc(meta_resp(bininds),edges);
%     [nvar] = histc(meta_var(bininds),edges);
%     sortRESP (:,ibin) = nresp;
%     sortVAR (:,ibin) = nvar;
% end