% sz = mdsize(A, dim)
% INPUTS:
%   A   = array
%   dim = dimensions
% OUTPUTS:
%   sz  = [size(A,dim(1)), size(A,dim(2)), ...]
function sz = mdsize(A, dims)
    for i = 1:length(dims)
        sz(i) = size(A, dims(i));
    end
end

