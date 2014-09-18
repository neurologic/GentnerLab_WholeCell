% operations for 2nd manuscript... intensity manipulation and response
% consistency

r=rigdef('mac')

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

%% get a list of all stimuli used across experiments (that have 5 rep or more)
trials = 5;
sig_wavname = [];
sig_wav = [];
tmpnames = [];
sigind = 1;
for iexpt=1:size(repexpts,1)
  
    thisexpt=repexpts{iexpt};
    load([r.Dir.Expt thisexpt])
    vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
    table=getClampTab(expt,{'clamp',0});
    keepsigs=reprequire(table,trials);
    thiscond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
    
    for isig = 1:size(thiscond,2)
        if ~isempty(regexp(thiscond(isig).wavnames,'d'))
        sig_wavname{sigind} = thiscond(isig).wavnames;
        tmpnames = strvcat(tmpnames ,thiscond(isig).wavnames);
        sigind = sigind + 1;
          end
    end
    
end

unq_wavname = unique(sig_wavname);
unq_used = [];
for iunq = 1:size(unq_wavname,2)
    numuse = 0;
   for iname = 1:size(sig_wavname,2)
      if ~isempty(regexp(sig_wavname{iname},unq_wavname{iunq}))
         
            numuse = numuse + 1;
          end
    
   end
   unq_used(iunq)=numuse;
end

for iexpt=1:size(repexpts,1)
  
    thisexpt=repexpts{iexpt};
    load([r.Dir.Expt thisexpt])
    vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
    table=getClampTab(expt,{'clamp',0});
    keepsigs=reprequire(table,trials);
    thiscond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
    
    for isig = 1:size(thiscond,2)
       
        sig_wavname{sigind} = thiscond(isig).wavnames;
        tmpnames = strvcat(tmpnames ,thiscond(isig).wavnames);
        sigind = sigind + 1;
    end
    
end

%% Vm distribution with changes in Vm
doplot = 1;

    stimind = 1;
for iexpt = 1:size(repexpts,1)
    thisexpt = repexpts{iexpt};
    load([r.Dir.Expt thisexpt]);
    
    vmexpt = filtesweeps(expt,0,'Vm',0);
    
    
    [sigon,sigoff]=GetSigTimes(vmexpt,vmexpt.stimcond,1);
    [dbstimcond,dblevels]=getDBstimcond(expt);

    vmedges = [-80:2:-30];
    for icond = 1:size(dbstimcond,1)
        thiscond = dbstimcond{icond};
        thisdb = dblevels{icond};
        vmcdf = [];
        vmhist = [];
        vmall = [];
        vmvar = [];
        vmsort = [];
        if doplot ==1 
        h_cdf = figure;
        hold on
        h_hist = figure;
        hold on
        h_var = figure;
        hold on
        h_sort = figure;
        hold on
        colgrad = colorGradient([0.5 0.5 0.5],[0,0,0],size(thisdb,2));
        end
        for idb = 1:size(thisdb,2);
            thiswav = thiscond(idb).wavnames;
            wavexpt = filtesweeps(vmexpt,0,'wavnames',thiswav);
            sigdata = medfilt1(wavexpt.wc.data*1000,200,[],2);
            sigdata = sigdata(:,sigon:sigoff);
            vmvar(idb) = mean(var(sigdata));
            
            sigdata = mean(sigdata);
            vmall(idb,:) = sigdata;
            
            vmsort(idb,:) = sort(sigdata);
             if doplot == 1
            figure(h_sort)
            line([1:size(sigdata,2)],vmsort(idb,:),'color',colgrad(idb,:),'LineWidth',4)
             end
             
            [x,p] = empcdf(sigdata);
            vmcdf(idb,:) = x;
            if doplot == 1
            figure(h_cdf)
            stairs(x,p,'color',colgrad(idb,:),'LineWidth',5)
            end
            
            n = histc(sigdata,vmedges);
            vmhist(idb,:) = n./size(sigdata,2);
            if doplot == 1
            figure(h_hist)
            stairs(vmedges,vmhist(idb,:),'color',colgrad(idb,:),'LineWidth',3)
            end

            
        end
        vmvar_expt{stimind} = vmvar;
        dblist_expt{stimind} = thisdb;
        
%         if doplot ==1
%         figure(h_var)
%         line(thisdb,vmvar./(max(vmvar)));
%         end
        
        stimind = stimind+1;
    end
    
end

