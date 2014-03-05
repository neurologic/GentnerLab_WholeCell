function [out_struct,hfig,hfig2,hfig3] = MetaResponseAnal_RsRin(allstepdata,expt)  
hfig2 = [];
stepdur = round(0.25/expt.wc.dt);
stepstart = 553;

stepdata = mean(allstepdata,1);
V_i = max(stepdata);
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

%need to implement fit options
options = fitoptions('exp2');
lowerTau = -1/expt.wc.dt;
upperTau = -1/0.05;
options.Lower = [0,lowerTau,0,lowerTau];
options.Upper = [Inf,upperTau,Inf,upperTau];
f=fit(xtimestep(1,Artifact_ind:100)',...
    stepdata(1,Artifact_ind:100)','exp2',options);
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
out_struct.TaoCell = rin_tau;
out_struct.Cm = (rin_tau / round(out_struct.Rin*1000000))*1000000000000; %seconds / ohms
% divide by 10^12 to return Capacitance in picoFarads

V_d = min(rin_fitline) + min(rs_fitline);
out_struct.V_f = V_i + V_d;

% subtract out electrode and non-voltage dependent component
tmp = stepdata - rs_fitline - rin_fitline;

%fit for slowly activating depol
tmp = tmp - max(tmp);
options.Lower = [-Inf,-Inf];
options.Upper = [0,0];
fd = fit(xtimestep', tmp' , 'exp1' ,options);
fitlineD = fd.a*exp(fd.b*xtimestep);
[rD rmse] = rsquare(tmp,fitlineD);
fitlineD = fitlineD - min(fitlineD);

hfig3 = figure;
hold on
line(xtimestep,tmp,'color','k','LineWidth',3)
line(xtimestep,fitlineD,'color','b','LineWidth',2)

% fit for slowly activating hyperpol
tmp = tmp - min(tmp);
options.Lower = [0,-Inf];
options.Upper = [Inf,0];
fh = fit(xtimestep', tmp' , 'exp1' ,options);
fitlineH = fh.a*exp(fh.b*xtimestep);
[rH rmse] = rsquare(tmp,fitlineH);
fitlineH = fitlineH - max(fitlineH);

%cant know conductance because can't know current
if rD >= rH 
    out_struct.Gv_dv = max(fitlineD);
    out_struct.TaoV = -1/fd.b;
    
    hfig2 = figure;
    hold on
    line(xtimestep,(stepdata- max(stepdata)+V_i),'color','k','LineWidth',2)
    line(xtimestep,rin_fitline+V_i,'color','g','LineWidth',2)
    line(xtimestep,rs_fitline+V_i,'color','r','LineWidth',2)
    line(xtimestep,fitlineD+V_i,'color','b','LineWidth',2)
end

if rH > rD
    out_struct.Gv_dv = min(fitlineD);
    out_struct.TaoV = -1/fh.b;
    
    hfig2 = figure;
    hold on
    line(xtimestep,(stepdata- max(stepdata)-V_i),'color','k','LineWidth',2)
    line(xtimestep,rin_fitline-V_i,'color','g','LineWidth',2)
    line(xtimestep,rs_fitline-V_i,'color','r','LineWidth',2)
    line(xtimestep,fitlineH-V_i,'color','b','LineWidth',2)
end



%% to plot when fooling around with the data bit by bit   
%     %%%%%%
% 
% figure;
% scatter(Rs,Rin,50,'k','fill');
% xlabel('Rs','FontSize',18),ylabel('Rin','FontSize',18)
% cf = fit(Rs',Rin','poly1');
% line(Rs, cf.p1.*Rs + cf.p2,'color','r')
% ci = confint(cf);
% title(['95% confidence interval on slope ('...
%     num2str(cf.p1) ') = [' num2str(ci(:,1)') ']'],'FontSize',18)
% 
% figure;
% scatter(Rs,Cm,50,'k','fill');
% xlabel('Rs','FontSize',18),ylabel('Cm','FontSize',18)
% cf = fit(Rs',Cm','poly1');
% line(Rs, cf.p1.*Rs + cf.p2,'color','r')
% ci = confint(cf);
% title(['95% confidence interval on slope ('...
%     num2str(cf.p1) ') = [' num2str(ci(:,1)') ']'],'FontSize',18)
% 
% figure;
% scatter(Rin,Cm,50,'k','fill');
% xlabel('Rin','FontSize',18),ylabel('Cm','FontSize',18)
% cf = fit(Rin',Cm','poly1');
% line(Rin, cf.p1.*Rin + cf.p2,'color','r')
% ci = confint(cf);
% title(['95% confidence interval on slope ('...
%     num2str(cf.p1) ') = [' num2str(ci(:,1)') ']'],'FontSize',18)
%
    


