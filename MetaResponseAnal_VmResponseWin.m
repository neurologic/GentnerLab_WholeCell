function MetaResponseAnal_VmResponseWin(expt,input_struct)
allfields = fieldnames(input_struct);
for ifield = 1:size(allfields,1)
   s = [allfields{ifield} ' = input_struct.' allfields{ifield} ';'];
   eval(s)
end

conf_p = 95;
windowsize = 500;
binsize = 10;
bin_p = 0.85;
ydatabound = [min(min(sigdata_filt(idb,basetimes(1):end))), max(max(sigdata_filt(idb,sigon:sigoff)))];
b=mean(sigdata_filt(idb,basetimes(1):basetimes(2)));
confint = getCDFconf (mean(sigdata_filt(:,basetimes(1):sigon)),conf_p);

[up_inds, low_inds] = WindowResponse(mean(sigdata_filt(:,sigon:sigoff)), confint, windowsize, binsize, bin_p);
up_win{sigind} = getWindowEdges (up_inds, 1, 1)+sigon;
low_win{sigind} = getWindowEdges (low_inds, 1, 1)+sigon;

for iresp = 1:size(up_win,2)
    spkinds = [];
    response_vm {sigind}(respind) = mean(mean(sigdata_filt(:,up_win(1,iresp):up_win(2,iresp))));
    for itrial=1:size(spiketimes,2)
        these_spk = spiketimes{itrial};
        spkinds = [spkinds, intersect(these_spk, [up_win(1,iresp):up_win(2,iresp)])];
    end
    response_spk {sigind}(respind) = (size(spkinds,2)/size(spiketimes,2))/(diff(up_win(:,iresp))*expt.wc.dt);
    respind = respind + 1;
end