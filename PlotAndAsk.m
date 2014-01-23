function outcell=PlotAndAsk(data,varargin)
 % data is vectors or matrix laid out along columns... 
 % so a vector of data is 1:alot... flip to plot
 
 % can only deal with asking for numerical data about the plot now
 
 % if asking more than one output use: cellfun(@eval,varargout) ...
 % when calling this function to actually assign your variables 
   %otherwise can just use eval(varargout)
 
 hfig=figure;
hold on
plot(data')

for iin=1:nargin-1
   tmpresp = input(varargin{iin});    
    outcell{iin}=[varargin{iin} '=' num2str(tmpresp) ';'];
end
nargin
close(hfig)