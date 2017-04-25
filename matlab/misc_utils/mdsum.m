% yNew = mdsum(y, dim)
%
% INPUTS:
%   y        = array
%   dim      = dimensions to sum along
% OUTPUTS:
%   ynew     = array with sums along dimensions in dim
function s = mdsum(y, dim)

s = y;
for d = dim
    s = sum(s,d);
end

end
