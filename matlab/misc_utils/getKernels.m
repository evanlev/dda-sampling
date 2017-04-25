% cell = getKernels(R, opt)
%
% R   = Reduction factor and period. 
%       For CAIPIRINHA, R = 8,10,12 => 15, 18, 28 possible kernels
% opt = 'all' => all kernels, 'caipirinha' => only caipirinha kernels
function cell = getKernels(R, opt)
    if nargin < 2
        opt = 'caipirinha';
    end
    opt = lower(opt);
    if strcmp(opt, 'caipirinha')
        cell = getCKernels(R);
    elseif strcmp(opt, 'all')
        cell = getPKernels([R R], R)
    else
        error('Option should be ''caipirinha'' or ''all'': %s', opt);
    end

end

function cell = getCKernels(R)
    cell = zeros(R,R,100);
     
    cell_tail = 1;
    for Ry = 1:R
    Rz = R / Ry;
    if (round(Rz) == Rz)
        for Delta = 0:Rz-1 %(R/Ry - 1)
            cell(1,1:Rz:end,cell_tail) = 1;
            last_row = cell(1,:,cell_tail);
            for s = Ry:Ry:R-1
                % Step down by Ry
                cell(1+s,:,cell_tail) = circshift(last_row, [0 Delta]);
                last_row = cell(1+s,:,cell_tail);
            end
            cell_tail = cell_tail + 1;
        end
    end
    end
    cell = cell(:,:,1:cell_tail-1); % Don't feel like thinking about it
end

% cells = getPKernels(cell_dims, k)
% Form all periodic sampling kernels, even the "bad" ones
%
% INPUTS:
%   cell_dims = [N M] 
%   k         = number of samples in a kernel
% OUTPUTS:
%   cells     = [N M P] array of kernels, cells(1,1,:) = 1
function sk_cells = getPKernels(cell_dims, k)
    RETURN_LIMIT = 1e7;
    assert(k >= 1);
    assert(k <= prod(cell_dims));
    if nchoosek(prod(cell_dims), k) > RETURN_LIMIT
        error('Cannot return more than %d kernels', RETURN_LIMIT);
    end
    
    cells = nchoosek(2:prod(cell_dims), k-1);
    inds = cells + repmat((0:size(cells,1)-1)*prod(cell_dims), [k-1 1])';
     
    sk_cells = zeros([cell_dims, size(cells,1)]);
    sk_cells(1,1,:) = 1;
    sk_cells(inds(:)) = 1;

end
