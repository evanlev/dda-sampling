% Cost reduction of adding
%
% c = getDeltaJ(pat, w, use_cpp)
% INPUTS:
%   pat    = sampling pat [ny nz]
%   w       = weighting function
%   use_cpp = 1/0 to use C++ implementation
% OUTPUTS:
%   c    = cost reduction from adding a sample 
%          [ny nz]
function deltaJ = getDeltaJ(pat, w, use_cpp)

fft_method = true; % Change to false to use the convolution method

if nargin < 3
    use_cpp = true;
end
for i = 1:3
    assert(size(pat,i) == size(w,i));
end
assert(size(w,3) == size(w,4));
assert(size(w,5) == 1 || size(w,5) == 2);

% Check inputs
if ndims(pat) > 3
    save pat_err_dims pat;
    error('pat must have at most 3 dims');
end

if use_cpp
    if ~fft_method % C++ version
        pat_file = [tempname '_egl'];
        w_file = [tempname '_egl'];
        out_file = [tempname, '_egl'];
        writeMDArray(int32(pat), pat_file);
        writeMDArray(w, w_file, 5);
        % Can add --debug <debug_level>
        args = [' --w ' w_file, ' --dJ ', out_file, ' --pat ', pat_file]; %, ' --w_type ', num2str(w_type), ' '];

        run_sysstr(['getDeltaJ  ', args]);

        deltaJ = readMDArray(out_file, 'double');

        % Clean up
        rm_dat(pat_file, w_file, out_file);
    else
        % FFT
        deltaJ = bart_0700('dda_getDeltaJ', w, pat);
        deltaJ = abs(deltaJ);
    end
else
    % MATLAB
    if fft_method && size(w,5) == 1
        % FFT method

        % TODO: partial Fourier
        nt = size(w,3);
        Fpat = fft(fft(pat, [], 1), [], 2);
        Fw = fft(fft(w, [], 1), [], 2);

        deltaJ = zeros(size(w,1), size(w,2), nt);
        for t1 = 1:nt
            deltaJ (:,:,t1) = w(1,1,t1,t1);
            for t2 = 1:nt
                deltaJ(:,:,t1) = deltaJ(:,:,t1) + ifft(ifft(Fpat(:,:,t2) .* (Fw(:,:,t1,t2) + circRev2d2(Fw(:,:,t2,t1))), [], 1), [], 2);
            end
        end
        deltaJ = abs(deltaJ);
    else
        if exist('getDeltaJ_mex') == 3
            %% MATLAB MEX
            assert(all(isreal(w(:))));
            deltaJ = getDeltaJ_mex(int32(pat), double(w));
        else
            fprintf('WARNING: do not have getDeltaJ_mex, reverting to matlab version...\n');
            %% MATLAB direct method
            if size(w,5) == 2
                deltaJ = getDeltaJ2(mdslice(w,5,1), pat, 1) ...
                       + getDeltaJ2(mdslice(w,5,2), pat, 2);
                deltaJ = deltaJ / 2;
            else
                deltaJ = getDeltaJ2(w, pat, 1);
            end
        end
    end
end % use MATLAB

deltaJ = real(deltaJ);

end % end getDeltaJ
 
% 2D circulant reversal of X
function Xr = circRev2d2(X)
    assert(ndims(X) <= 2);
    Xr = zeros(size(X));
    Xr(2:end,1) = flipdim(X(2:end,1),1);
    Xr(1,2:end) = flipdim(X(1,2:end),2);
    Xr(2:end,2:end) = flipdim(flipdim(X(2:end,2:end),1),2);
    Xr(1,1) = X(1,1);
end

function run_sysstr(comm)
    fprintf('executing: %s\n', comm);
    system(comm);
end

function rm_dat(varargin)
    for i = 1:length(varargin)
        delete([varargin{i} '.dat']);
        delete([varargin{i} '.hdr']);
    end
end

function deltaJ = getDeltaJ2(w, pat, w_type)
    % Mask dimensions
    [ny nz nt] = size(pat);
    assert(size(w,5) == 1);

    % Compute deltaJ sample add cost
    fprintf('Computing add cost\n');
    deltaJ = zeros(ny,nz,nt);
    for t = 1:nt
        deltaJ(:,:,t) = w(1,1,t,t);
    end
    [sy, sz, st] = ind2sub(size(pat), find(pat(:)~=0));
    sy = sy - 1;
    sz = sz - 1; 
    st = st - 1; % C indexing
    Ns = sum(pat(:));
    for i = 1:Ns
        for t_new  = 0:nt-1
        for ky_new = 0:ny-1
        for kz_new = 0:nz-1
            % First part
            if w_type == 1
                dy = mod(ky_new - sy(i), ny)+1;
                dz = mod(kz_new - sz(i), nz)+1;
            else
                dy = mod(sy(i) + ky_new, ny)+1;
                dz = mod(sz(i) + kz_new, nz)+1;
            end
            %fprintf('[%d %d] - [%d %d] = [%d %d %d] ', sy(i), sz(i), ky_new, kz_new, dy, dz, t_new);
            deltaJ(ky_new+1,kz_new+1,t_new+1) = deltaJ(ky_new+1,kz_new+1,t_new+1) ...
                        + w(dy,dz,t_new+1,st(i)+1);
            if w_type == 1
                dy = mod(sy(i) - ky_new, ny)+1;
                dz = mod(sz(i) - kz_new, nz)+1;
            end
            deltaJ(ky_new+1,kz_new+1,t_new+1) = deltaJ(ky_new+1,kz_new+1,t_new+1) ...
                        + w(dy,dz,st(i)+1,t_new+1);
        end
        end
        end
        fprintf('\r%d%%', round(i/Ns * 100));
    end
end


