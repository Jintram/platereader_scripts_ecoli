function dispAndWrite(fid,text)%,varargin)
%
%

%{
if nargin>2
    for i=1:2:(numel(varargin)-1)
        internalSettings.(varargin{i}) = varargin{i+1};
    end
else
    internalSettings=struct;
end

if ~isfield(internalSettings, 'noDisp');
    disp(text); 
elseif ~internalSettings.noDisp
    disp(text); 
end
%}

disp(text); 
fprintf(fid, '%s\n', text); 