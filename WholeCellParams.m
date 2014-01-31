function [Rin,Rs]=WholeCellParams(expt,isIC)
% provide an expt that is filtered to only one holding potential / clamp
%0.014 seems like the time to wait for the RC of the Rs transient before
%measuring Rin
if isIC==1
    Rs=[];
    Rin=[];
    stepampcommand=-75; % in picoAmps
     stepdata=mean(expt.wc.data(:,1:expt.analysis.params.baselinewin(2)),1);
    baselineval=median(stepdata(1,expt.analysis.params.baselinewin(1):expt.analysis.params.baselinewin(2)));
    stepval=median(stepdata(1,expt.analysis.params.steptime(1)+round((0.014/expt.wc.dt)):expt.analysis.params.steptime(2)));
    dv=diff([baselineval,stepval]);
    Rin=round((dv/stepampcommand)*10^6); %     values in current clamp are recorded in Volts
    
end

if isIC==0;
    stepampcommand=-10;
     stepdata=mean(expt.wc.data(:,1:expt.analysis.params.baselinewin(2)),1);
    baselineval=median(stepdata(1,expt.analysis.params.baselinewin(1):expt.analysis.params.baselinewin(2)));
    stepval=median(stepdata(1,expt.analysis.params.steptime(1)+round((0.014/expt.wc.dt)):expt.analysis.params.steptime(2)));
    dA=diff([baselineval,stepval]);
    Rin=round((stepampcommand/dA)*10^3);
    
    transval=min(stepdata(1,1:expt.analysis.params.steptime(1)+round((0.014/expt.wc.dt))));
     dA=diff([baselineval,transval]);
    Rs=round((stepampcommand/dA)*10^3);
end

