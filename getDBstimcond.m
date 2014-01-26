function [dbstimcond,dblevels]=getDBstimcond(expt)

% filter stimcond for DB signals... 
% get out how many unique there are
% how many DB and what DB of each
stimcond=expt.stimcond;
thesesigs=[];
dbwavnames=[];
dblevels=[];
stimind=1;
for istim=1:size(expt.stimcond,2)
    dbind=regexp(stimcond(istim).wavnames,'d');
    if ~isempty(dbind)
        thesesigs(istim)=1;
        dbwavnames{stimind}=stimcond(istim).wavnames(1:dbind-1);
        dblevelsall(stimind)=str2num(stimcond(istim).wavnames(end-1:end));
        stimind=stimind+1;
    else thesesigs(istim)=0;
    end
    
end
dbsigs=find(thesesigs);
stimcond=expt.stimcond(dbsigs);

dbstims=unique(dbwavnames);
dbstimcond=[];
if isempty(dbstims)
    return
end
for istim=1:size(dbstims,2)
    for iall=1:size(dbwavnames,2)
        %get indices for all instinces of that wav
        iswav=regexp(dbwavnames{iall},dbstims{istim});
        if ~isempty(iswav)
            dbinds(iall)= 1;
        else dbinds(iall)= 0;
        end
    end
    tmpcond=stimcond(find(dbinds));
    dblevels{istim}=dblevelsall(find(dbinds));
    dbstimcond{istim}=tmpcond;
end