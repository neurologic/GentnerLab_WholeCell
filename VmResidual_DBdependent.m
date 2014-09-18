%for each experiment, get the unique stimuli
%for each stimulus see if there were multiple dB played
%for each dB get skewness of Vm residuals
stimind = 1;
for iexpt=1:size(repexpts,1)
    thisexpt=repexpts{iexpt};
    load([r.Dir.Expt thisexpt])
    vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
    table=getClampTab(expt,{'clamp',0});
    thiscond=getsubstimcond(expt.stimcond,table.sigsplayed);
    trials = 5;
    keepsigs=reprequire(table,trials);
    allsig=table.sigsplayed;
    repsigs=allsig(keepsigs);
    repstimcond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
    
    [dbstimcond,dblevels]=getDBstimcond(expt);
   
    for istim = 1:size(dbstimcond)
         dbaxis{stimind} = dblevels{istim};
         dbcond = [];
         dbnames = [];
        for idb = 1:size(dbstimcond{istim},2)
            dbcond = dbstimcond{istim};
            dbnames{idb} = dbcond(idb).wavnames;
        end
        stimcond=getsubstimcond(repstimcond,dbnames');
        skew_db = [];
        var_db = [];
        for idb = 1:size(stimcond,2)
            [sigon,sigoff]=GetSigTimes(vmexpt,stimcond,1);
            basetimes = vmexpt.analysis.params.baselinewin;
            sigexpt = filtesweeps(vmexpt,0,'wavnames',stimcond(idb).wavnames);
            sigdata = sigexpt.wc.data*1000;
            sigdata = medfilt1(sigdata,200,[],2);
            sigdata = sigdata(:,sigon:sigoff);
            var_db(idb) = mean(var(sigdata,1),2);
            residat = (sigdata-repmat(mean(sigdata,1),size(sigdata,1),1));
            
            skew_db(idb) = skewness(residat(:));
        end
        skewaxis{stimind} = skew_db;
        varaxis{stimind} = var_db;
        stimind = stimind+1;
    end
end

%%
figure;
hold on
for istim = [1:9,11:20]
    scatter(dbaxis{istim},skewaxis{istim})
end

figure;
hold on
for istim = [1:9,11:20]
    scatter(dbaxis{istim},varaxis{istim}/max(varaxis{istim}))
end