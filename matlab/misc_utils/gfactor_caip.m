% g = gfactor_caip(sns_maps, kernel)
%
% g-factor for CAIPIRINHA
%
% Evan Levine, 4/2/2016
%
% INPUTS:
%   sns_maps = sensitivity maps
%   kernel   = period of the sampling pattern
% 
% References:
%   Breuer, Felix A., et al. "Controlled aliasing in volumetric parallel 
%   imaging (2D CAIPIRINHA)." Magnetic resonance in medicine 55.3 (2006): 
%   549-556.
function gmap = gfactor_caip(sns_maps, kernel, psi)

    tic;

    % Check input
    assert(ndims(sns_maps) == 4); % TODO multiple maps 
    assert(ndims(kernel) == 2);

    % Size
    [ny0 nz0] = size(kernel);
    [nx ny nz ncoils] = size(sns_maps);
    
    if nargin < 3
        psi = eye(ncoils);
    end

    % Get aliases
    aliases = fftn(kernel);
    aliases(abs(aliases) < 0.01) = 0;
    aliases = find(aliases(:));
    [ay, az] = ind2sub([ny0 nz0], aliases);
    gmap = zeros(nx,ny,nz);

    % Store psi^-1
    psi_inv = inv(psi);

    % g-factor
    gmap = zeros(nx,ny,nz);
    for x = 1:nx
    for y = 1:ny
    for z = 1:nz
        % S matrix and aliases
        S = zeros(ncoils, length(aliases));
        ya = zeros(length(aliases),1);
        za = zeros(length(aliases),1);
        for i = 1:length(aliases)
            % Voxels aliasing into this one
            ya(i) = mod(y + round(ny*ay(i)/ny0)-1, ny)+1;
            za(i) = mod(z + round(nz*az(i)/nz0)-1, nz)+1;
            S(:,i) = squeeze(sns_maps(x,ya(i),za(i),:));
        end

        % Remove columns corresponding to voxels outside of the support
        support = (sum(S ~= 0, 1) ~= 0); 
        S = S(:,support);
        ya = ya(support);
        za = za(support);
        aliases_yz = aliases(support);

        % g-factor from original SENSE paper
        gsq = zeros(length(aliases_yz),1);
        if ~isempty(S)
            StPS = S'*psi_inv*S;
            g = sqrt(diag(inv(StPS)) .* diag(StPS));
            for i = 1:length(gsq)
                gmap(x,ya(i),za(i)) = g(i);
            end
        end
    end % end z loop
    end % end y loop
    end % end x loop
    gmap = real(gmap);

    fprintf('g-factor for caip done in %f sec\n', toc);
end

