function substimcond=getsubstimcond(stimcond,sublist)

% sublist=unique(Vmexpt.sweeps.wavnames);
metalist=[];
for iwav=1:size(sublist,1)
    metalist=[metalist, sublist(iwav)];
end
stiminds=[];
for istimcond=1:size(stimcond,2)
    metainds=[];
    for imeta=1:max(size(metalist))
        if ~isempty(strfind(metalist{imeta},stimcond(istimcond).wavnames))
            if size(metalist{imeta},2)==size(stimcond(istimcond).wavnames,2)
                metainds(imeta)= 1;
            end
        end
        if isempty(strfind(metalist{imeta},stimcond(istimcond).wavnames))
            metainds(imeta)=0;
        end
    end
    
    if ~isempty(find(metainds))
        stiminds(istimcond)=1;
    end
   if isempty(find(metainds))
       stiminds(istimcond)=0;
    end
end
keepstims=find(stiminds);
substimcond=stimcond(keepstims);

