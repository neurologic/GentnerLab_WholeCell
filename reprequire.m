function keepsigs=reprequire(table,numreps)
%works off of a very specific table struct created using QueryExpt
% the numtrials field of this struct (which is saved in an expt) is queried
for isig=1:max(size(table.numtrials))
   if table.numtrials(isig) < numreps;
       excluded(isig)=0;
   end
   if table.numtrials(isig) >= numreps;
       excluded(isig)=1;
   end
end
keepsigs=find(excluded);