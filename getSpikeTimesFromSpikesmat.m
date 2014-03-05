function spikeInds = getSpikeTimesFromSpikesmat(spikesmat,highpassdata)
for itrial=1:size(spikesmat,1)
    trial_data = highpassdata(itrial,:);
    spk_inds = find(spikesmat(itrial,:));
    tmp_inds = zeros(1,size(spikesmat,2));
    if ~isempty(spk_inds)
        for ispike = 1:size(spk_inds)
            this_win = [spk_inds(ispike),spk_inds(ispike)+20];
            tmp_data = trial_data(1,this_win(1):this_win(2));
            peak = min(find(tmp_data == max(tmp_data)));
            peak_inds(ispike) = spk_inds(ispike)+peak;
        end
        tmp_inds(1,peak_inds) = 1;
    end
    spikeInds(itrial,:) = tmp_inds;
end