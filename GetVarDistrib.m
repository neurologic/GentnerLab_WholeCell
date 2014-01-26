function [sigdata,vardata] = GetVarDistrib(expt,stimcond)

for isig=1:size(stimcond,2)
    sigexpt=filtesweeps(expt,0,'wavnames',stimcond(isig).wavnames);
    filtdata=medfilt1(sigexpt.wc.data,200,[],2);
    
    sigdata(isig,:)=mean(filtdata);
    vardata(isig,:)=var(filtdata);
    varallowt=round(0.1/expt.wc.dt); %use this filter to allow variance in real time to go above signif a little
    
    vardataf=vardata(isig,:); % dont filter so much anymore... just start looking for offset after ~
    
end

