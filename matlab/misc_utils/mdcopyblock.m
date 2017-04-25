% yNew = mdcopyblock(y, dim, block)
% INPUTS:
%   y     = array
%   dim   = dimensions to extract block
%   block = range to extract (e.g. 2:4), if length(dim) > 1,
%           provide a cell array of ranges, otherwise no need
% OUTPUTS:
%   yNew = array(:,:,block1,:,:,block2,:,:)
%                     ^           ^
%                   dim(1)      dim(2)  ....
%
% See also: mdslice, mdrepmat, mdcrop, mdzpad
function yNew = mdcopyblock(y, dim, block)
    if ~(length(dim) == length(block) || ...
            (length(dim) == 1 && ~iscell(block)))
        dim
        block
        error('incompatible dim/block above');
    end
    if iscell(block)
        for i = 1:length(dim)
            if ~(size(y,dim(i)) >= length(block{i}))
                dim(i)
                size(y,dim(i))
                block{i}
                i
                length(block{i})
                error('Bad dim/size/block/i above');
            end
        end
    end
    if length(dim) == 1 && iscell(block)
        block = block{1};
    elseif length(dim) > 1
        for i = 1:length(dim)
            y = mdcopyblock(y, dim(i), block{i});
        end
        yNew = y;
        return;
    end
    M = size(y, dim);
    prod_other_dims = numel(y) / M;

    perm = 1:ndims(y);
    perm(1) = dim;
    perm(dim) = 1;
    
    yNew = permute(y, perm);
    perm_size = size(yNew);
    yNew = reshape(yNew, [M, prod_other_dims]);
    yNew = yNew(block,:);
    perm_size(1) = size(yNew, 1);
    yNew = reshape(yNew, perm_size);
    yNew = ipermute(yNew, perm);
end


