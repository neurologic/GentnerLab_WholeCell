function [out_struct,hfig] = MetaResponseAnal_RsRin(expt,input_struct)
allfields = fieldnames(input_struct);
for ifield = 1:size(allfields,1)
    s = [allfields{ifield} ' = input_struct.' allfields{ifield} ';'];
    eval(s)
end

 
stepdur = round(0.25/expt.wc.dt);
 stepstart = 553;
allstepdata = expt.wc.data(:,stepstart:stepdur)*1000;
stepstd = std(allstepdata');

approxDV = min(min(expt.wc.data(:,1:stepstart)))*1000 - min(min(allstepdata));
%only use trials to calculate parameters that have standard deviation less
%than half the approximated voltage drop during step (not too much spont
%activity riding on step)
keepinds = find(stepstd < approxDV/2);

stepdata = mean(allstepdata(keepinds,:),1);
stepdata=(stepdata- max(stepdata));
rin_x = [0:size(stepdata,2)-1];

hfig = figure; hold on
xtimestep = rin_x*expt.wc.dt;
line(xtimestep,stepdata,'color','k','LineWidth',2)

dV = min(stepdata) 
stepdata = stepdata - dV;

if ~isempty(find(diff(stepdata(1,1:10))>0))
    Artifact_ind = 4;
elseif isempty(find(diff(stepdata(1,1:10))>0))
    Artifact_ind = 1;
end


f=fit(xtimestep(1,Artifact_ind:100)',stepdata(1,Artifact_ind:100)','exp2');
alltau= [f.b,f.d];
alloff = [f.a,f.c];
rs_t = min(alltau);
rs_off = alloff(find(alltau == rs_t));
rin_t = max(alltau);
rin_off = alloff(find(alltau == rin_t));
rs_fitline = rs_off*exp(rs_t*xtimestep);
rs_fitline = rs_fitline - max(rs_fitline);
line(xtimestep,rs_fitline,'color','r');
rin_fitline = rin_off*exp(rin_t*xtimestep);
rin_fitline = rin_fitline - max(rin_fitline);
line(xtimestep,rin_fitline,'color','g');
line(xtimestep,(rin_fitline+rs_fitline),'color','b')

stepamp=-0.075;
out_struct.Rin  = min(rin_fitline)/stepamp;
out_struct.Rs = min(rs_fitline)/stepamp;
rin_tau = -1/rin_t;  
out_struct.Cm = (rin_tau / round(out_struct.Rin*1000000))*1000000000000; %seconds / ohms
    % divide by 10^12 to return Capacitance in picoFarads
