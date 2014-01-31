function MetaResponseAnal_VmResiduals(expt,input_struct)
allfields = fieldnames(input_struct);
for ifield = 1:size(allfields,1)
   s = [allfields{ifield} ' = input_struct.' allfields{ifield} ';'];
   eval(s)
end

resid_data = sigdata_filt - repmat(mean(sigdata_filt),size(sigdata_filt,1),1);
resid_edges = [-40:2:40];
tmp_stim = resid_data(:,sigon:sigoff);
tmp_stim = reshape(tmp_stim',1,...
    size(tmp_stim,1)*size(tmp_stim,2));
tmp_base = resid_data(:,basetimes(1):sigon);
tmp_base = reshape(tmp_base',1,...
    size(tmp_base,1)*size(tmp_base,2));
nresid_stim(stimind,:) = histc(tmp_stim,resid_edges)./size(tmp_stim,2);
nresid_base(stimind,:) = histc(tmp_base,resid_edges)./size(tmp_base,2);

skew_base (stimind) = skewness(tmp_base);
skew_stim (stimind) = skewness(tmp_stim);

if plotResidVm ==1
    figure;
    hold on
    stairs(resid_edges,nresid_base(stimind,:),'color','b','LineWidth',3)
    stairs(resid_edges,nresid_stim(stimind,:),'color','r','LineWidth',3)
    scatter(mean(tmp_stim),max([nresid_stim(stimind,:),nresid_base(stimind,:)])+0.05,100,'r','v')
    scatter(mean(tmp_base),max([nresid_stim(stimind,:),nresid_base(stimind,:)])+0.05,100,'b','v')
    ylims = get(gca,'YLim');
    text(resid_edges(10),ylims(2)-0.05,{['base skew = ' num2str(skew_base(stimind))];...
        ['stim skew = ' num2str(skew_stim(stimind))]},'HorizontalAlignment','center');
    %using residuals is a really easy way to "highpass"
    title([expt.name ';  stimulus# ' num2str(istimcond)...
        ';  distribution Vm residuals;  stim = ' ...
        num2str(thisdb(dbind)) 'dBSPL'],'Interpreter','none');
end