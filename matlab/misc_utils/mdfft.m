% yNew = mdfft(y, dim, unitary, inv, [centered])
%
% INPUTS:
%   y        = array
%   dim      = dimensions to apply FFT (e.g. [1 2 4 7])
%   unitary  = normalize
%   inv      = apply inverse FFT
%   centered = centered (default=true)
% OUTPUTS:
%   ynew    = Fourier-transformed array
function y = mdfft(y, dim, unitary, inv, centered)
    if nargin < 5
        centered = true;
    end
    if length(dim) > 1
        % Don't use fftn for simplicity
        for i = 1:length(dim)
            if size(y, dim(i)) > 1
                y = mdfft(y, dim(i), unitary, inv, centered);
            end
        end
        return;
    else
        % Apply centered FFT 
        if inv
            if centered
                y = ifftshift(ifft(ifftshift(y, dim), [], dim), dim);
            else
                y = ifft(y, [], dim);
            end
        else
            if centered
                %fprintf('fft along %d\n', dim);
                y = fftshift(fft(fftshift(y, dim), [], dim), dim);
            else
                y = fft(y, [], dim);
            end
        end 
        % Normalize
        if unitary
            if inv
                %fprintf('x %f\n', prod(mdsize(y,dim)));
                y = y * sqrt(prod(mdsize(y, dim)));
            else
                %fprintf('Divide by %f\n', prod(mdsize(y,dim)));
                y = y / sqrt(prod(mdsize(y, dim)));
            end
        end
    end
end

