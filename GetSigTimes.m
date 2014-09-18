function [sigon,sigoff]=GetSigTimes(expt,stimcond,isig)
sigon=round(expt.analysis.params.waveonset_time/expt.wc.dt);

%for sigoff need to check if it it is a short stim that was padded...
tmpy = stimcond(isig).wavs;
%check if last 50 points is all zeros
sigoff=round(sigon+((length(stimcond(isig).wavs)/44100)/expt.wc.dt));
