% yNew = mdcrop(y, dim, s)
%
% INPUTS:
%   y   = array of any size
%   dim = dimension to crop
%   s   = size of 
% OUTPUTS:
%   yNew = array y cropped around its center
%
% See also: mdslice, mdrepmat, mdcopyblock, mdzpad
function yNew = mdcrop(y, dim, s)
    for i = 1:length(dim)
        assert(s(i) <= size(y,dim(i)));
    end
    if length(dim) > 1
        for i = 1:length(dim)
            y = mdcrop(y, dim(i), s(i));
        end
        yNew = y;
        return;
    end

    M = size(y,dim);
	range = floor(M/2)+1+ceil(-s/2) : floor(M/2)+ceil(s/2);
    assert(length(range) == s);
    assert(all(range) >= 1);
    assert(all(range) <= M);
    yNew = mdcopyblock(y, dim, range);
end

