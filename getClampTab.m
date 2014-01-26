function table=getClampTab(expt,varargin);
%varargin is in the format: {field, type}
fields=fieldnames(expt.table);
filtinds=[];
for iin=1:size(varargin,2)
    thisfield=varargin{iin}{1};
    thisval=varargin{iin}{2};
    
    for ifield=1:size(fields ,1)
        if ~isempty(regexp(fields{ifield},thisfield))
            if size(thisfield,2)==size(fields{ifield},2)
                for ifilt=1:size(expt.table,2)
                    s=['tmp=expt.table(' num2str(ifilt) ').' thisfield ';'];
                    eval(s)
                    if ~ischar(thisval)
                        for ival=1:max(size(tmp))
                            if ~isempty(find(tmp(ival)==thisval))
                                filtinds=[filtinds, ifilt];
                               
                            end
                        end
                    end
                    if ischar(thisval)
                        for ival=1:max(size(tmp))
                            if ~isempty(regexp(tmp{ival},thisval))
                                filtinds=[filtinds, ifilt];
                                
                            end
                        end
                    end
                end
            end
        end
    end
end




table=expt.table(unique(filtinds));
% for