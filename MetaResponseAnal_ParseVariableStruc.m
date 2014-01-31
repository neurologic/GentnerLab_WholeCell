function varargout = MetaResponseAnal_ParseVariableStruc(input_struct)

allfields = fieldnames(input_struct);

for ifield = 1:size(allfields,1)
   s = [allfields{ifield} ' = input_struct.' allfields{ifield} ';'];
%    eval(s)
   varargout{ifield} = s;
end