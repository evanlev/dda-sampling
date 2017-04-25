% res = mdslice(array, dim, pos)
%
% Slice an array 
% INPUTS:
%   y     = array
%   dim   = [D 1] array of dimensions to slice
%   pos   = [D 1] array of slice positions
% OUTPUTS:
%   res = array[:,:,:,pos,:,:,:]
%                      ^
%                dimension dim
%
% See also: mdslice, mdrepmat, mdcrop, mdzpad
function array_perm = mdslice(array, dim, pos)
    % Check input
    assert(all(pos >= 1));
    assert(all(dim(pos > 1) <= ndims(array)));
    for i = 1:length(dim)
        assert(pos(i) <= size(array,dim(i)));
    end

    % slice dim-by-dim
    array_perm = array;
    for i = 1:length(dim)
        array_perm = mdslice1(array_perm, dim(i), pos(i));
    end
end

% mdslice with 1 dim
function array_perm = mdslice1(array, dim, pos)
    assert(length(dim) == 1 && length(pos) == 1);

    prod_other_dims = numel(array) / size(array, dim);
 
    perm = 1:max(dim, ndims(array));
    perm(1) = dim;
    perm(dim) = 1;
    
    array_size = size(array);
    array_perm = permute(array, perm);
    array_size_perm = size(array_perm);
    array_perm = reshape(array_perm, [size(array,dim), prod_other_dims]);

    array_perm = array_perm(pos,:);
    array_size_perm(1) = 1; 
    
    array_perm = ipermute(reshape(array_perm, array_size_perm), perm);

end
