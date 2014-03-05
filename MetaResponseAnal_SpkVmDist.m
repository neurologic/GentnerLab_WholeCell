function out_struct = MetaResponseAnal_SpkVmDist(expt,input_struct)
allfields = fieldnames(input_struct);
for ifield = 1:size(allfields,1)
    s = [allfields{ifield} ' = input_struct.' allfields{ifield} ';'];
    eval(s)
end

[spikesmat]=getspikesmat(highpassdata,spk_thresh,expt.wc.dt);
spikeInds = getSpikeTimesFromSpikesmat(spikesmat,highpassdata);

out_struct.VmBins = histc(sigdata_filt(:),VmEdges)./size(sigdata_filt(:),1);
out_struct.numVm = size(sigdata_filt(:),1);
spk_ind = 1;
clipped_vm = [];
for itrial = 1:size(spikeInds,1)
    these_spk = find(spikeInds(itrial,:));
    if ~isempty(these_spk)
        for ispike = 1:size(these_spk)
            clipped_vm(1,spk_ind) = sigdata_filt(itrial,these_spk(ispike));
            spk_ind = spk_ind + 1;
        end
    end
end
out_struct.SpkVmBins = histc(clipped_vm,VmEdges)./size(clipped_vm,2);
out_struct.numSpks = size(clipped_vm,2);

%do actual laborious p(spike|Vm)
binlen = 0.02; %seconds
binsize = round(binlen / dt);
sigdata_filt = sigdata_filt(:,1:size(spikeInds,2));
long_spk = spikeInds(:);
long_vm = sigdata_filt(:);
numbins = floor(size(long_spk,1)/binsize)
bins = [1:binsize:binsize*numbins];
numspk = zeros(1,numbins - 1);
for ibin = 1:numbins-1
   bin_vm(ibin) = mean(long_vm(bins(ibin):bins(ibin+1),1));
   is_spk = find(long_spk(bins(ibin):bins(ibin+1),1) == 1);
   if ~isempty(is_spk)
        numspk(1,ibin) = max(size(is_spk));
   end
end
[cnt_vm,ind] = histc(bin_vm,VmEdges);

cnt_spk = zeros(size(cnt_vm,1),size(cnt_vm,2));
for iedge = 1:size(VmEdges,2)
    these_ind = find(ind == iedge);
    these_spk = find(numspk(these_ind));
    if ~isempty(these_spk)
    cnt_spk(1,iedge) = max(size(these_spk));
    end
end

pVm = cnt_vm / size(bin_vm,2);
% pSpk = pVm .* (cnt_spk ./ cnt_vm);
out_struct.pSpk = (cnt_spk ./ cnt_vm);

    


