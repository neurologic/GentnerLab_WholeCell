function [isclamp,tableID]=findclamp(expt,clamptype)
% works on the "table" field of expt created using QueryExpt
% returns 1 if expt has that clamptype
% returnds 0 if the expt does not have that clamptype

isclamp=0;
tableID=0;
clampind=1;
for iclamp=1:size(expt.table,2)
    if regexp(expt.table(iclamp).clamptype,clamptype)
        isclamp=1;
        tableID(clampind)=iclamp;
        clampind=clampind+1;
    end
end