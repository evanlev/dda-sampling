%  p = get_dd(mask, [opt])
% INPUTS:
%   mask   = sampling pattern, ny x nz x nt
%   opt    = 'Real'    => use the real signal model
%            'Complex' => use the complex signal model (default)
% OUTPUTS:
%   p      = differential distribution
function p = get_dd(mask, opt, use_bart)

    % Check input
    assert(ndims(mask) <= 3);

    % Real / complex signal model
    if nargin < 2
        opt = 'COMPLEX';
    end
    if nargin < 3
        use_bart = true; % For testing
    end

    opt = upper(opt);
    assert(strcmp(opt, 'REAL') || strcmp(opt, 'COMPLEX'));
    
    if use_bart
        % Use BART O(samples ^ 2)
        p = bart_0700('dda_getp', mask);
    else
        % Use an FFT method, usually faster if the pattern is not sparse
        PSF = ifft(ifft(mask, [], 1), [], 2);
        p = get_dd2(PSF, opt);
    end
end

function p = get_dd2(PSF, opt)
    % Compute the first part
    [ny nz nt] = size(PSF);
    N = ny*nz;
    p = zeros(ny,nz,nt,nt);
    for t1 = 1:size(PSF,3)
    for t2 = 1:size(PSF,3)
        p(:,:,t1,t2) = N*fft(fft(PSF(:,:,t1) .* conj(PSF(:,:,t2)), [], 1), [], 2);
    end
    end

    % For real signal model, add the second part
    if strcmp(opt, 'REAL')
        p2 = zeros(ny,nz,nt,nt);

        fft2d = @(X) fft(fft(X, [], 1), [], 2);
        for t1 = 1:size(PSF,3)
        for t2 = 1:size(PSF,3)
            p2(:,:,t1,t2) = N * fft2d(PSF(:,:,t1).^2); % .* conj(flip2d(PSF(:,:,t2))));
        end
        end
        p = cat(4, p, p2);
    end
    
    % Imaginary part is due to finite precision
    assert(norm(vec(imag(p))) / norm(vec(real(p))) < 1e-7);
    p = real(p);
end


