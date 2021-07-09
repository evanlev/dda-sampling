%  W = buildW(sns_maps)
%
% INPUTS:
%   sns_maps = [Nx Ny Nz C M T L] array
% OUTPUTS:
%   w = [Ny Nz T T] array
% 
% Build w weighting such that <w, p> = sum_k \lambda_k(E'E)^2
function weighting = buildW(sns_maps, use_bart)
    if nargin < 2
        use_bart = true;
    end
    assert(ndims(sns_maps) < 9);

    fprintf('Building weighting...\n');
    if use_bart
        try
            weighting = buildW2(sns_maps, true);
        catch e
            disp('WARNING: Did not find BART, reverting to MATLAB version...\n');
            weighting = buildW2(sns_maps, false);
        end
    else
        weighting = buildW2(sns_maps, false);
    end
end

function weighting = buildW2(sns_maps, use_bart)

    if use_bart
        weighting = bart_0700('dda_getw', sns_maps);
    else
        TEMAP_DIM = 7;
        TE_DIM    = 6; % time (t variable)
        MAPS_DIM  = 5;
        COIL_DIM  = 4;
        READ_DIM  = 1;
    
        sns_dims = mdsize(sns_maps, 1:7);
        Nw = size(sns_maps,8);

        n_loop = prod(sns_dims([TE_DIM, TE_DIM, COIL_DIM, COIL_DIM, READ_DIM]))*Nw;
        loop_count = 0;

        w_dims = [sns_dims(2:3), sns_dims(TE_DIM), sns_dims(TE_DIM), Nw];
        weighting = zeros(w_dims);
        for ri = 1:Nw
        for t1 = 1:sns_dims(TE_DIM)
        for t2 = 1:sns_dims(TE_DIM)
        for c1 = 1:sns_dims(COIL_DIM)
        for c2 = 1:sns_dims(COIL_DIM)
        for x  = 1:sns_dims(READ_DIM)
            if ri == 1
                cj = @(X) conj(X);
            else
                cj = @(X) X;
            end
            LHS1 = abs(fftn(sum(sum(cj(sns_maps(x,:,:,c2,:,t2,:)).*...
                                       sns_maps(x,:,:,c1,:,t1,:), ...
                                       MAPS_DIM),TEMAP_DIM))).^2;

            weighting(:,:,t1,t2,ri) = weighting(:,:,t1,t2,ri) + shiftdim(LHS1, 1);

            loop_count = loop_count + 1;
            if( n_loop > 10 && mod(loop_count , floor(n_loop / 10)) == 0 )
                fprintf('\r%02d%%', round(100*loop_count / n_loop));
            end
        end % x
        end % c2
        end % c1
        end % t2
        end % t1
        end % ri
        weighting = 1/(Nw * prod(sns_dims(1:3))^2)*weighting;
    end

    % Should be real
    assert(norm(imag(weighting(:)))/norm(real(weighting(:))) < 1e-7);
    weighting = real(weighting);

end

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


