function [sigon,sigoff]=GetSigTimes(expt,stimcond,isig)
sigon=round(expt.analysis.params.waveonset_time/expt.wc.dt);
sigoff=round(sigon+((length(stimcond(isig).wavs)/44100)/expt.wc.dt));
