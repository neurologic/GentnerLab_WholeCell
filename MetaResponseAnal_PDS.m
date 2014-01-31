function MetaResponseAnal_PDS(expt, input_struct)

allfields = fieldnames(input_struct);
for ifield = 1:size(allfields,1)
   s = [allfields{ifield} ' = input_struct.' allfields{ifield} ';'];
   eval(s)
end

Fs = 1/expt.wc.dt;
s = sigdata(:,sigon:sigoff)';
L = size(s,1);
p = abs(fft(s))/(L/2); %% absolute value of the fft
p = p(1:L/2).^2 %% take the power of positve freq. half
freq = [0:L/2-1]/expt.wc.dt; %% find the corresponding frequency in Hz
semilogy(freq,p); %% plot on semilog scale
