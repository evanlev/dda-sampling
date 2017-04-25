% l1eigs_analytical = get_l1eigs(sns_maps, Nsamps)
%   INPUTS:
%       sns_maps          = sensitivity maps [nx ny nz ncoils nmaps nt ntmap]
%       Nsamps            = Number of samples for each echo [nt]
%   OUTPUTS:
%       l1eigs_analytical = scalar, sum of eigenvalues of E^HE, where E is the 
%                           SENSE encoding matrix
function l1eigs_analytical = get_l1eigs(sns_maps, Nsamps)
    TIME_DIM = 6;
    TIMEMAP_DIM = TIME_DIM + 1;
    l1eigs_analytical = 0;
    % Nx * Ny *nz
    N = numel(sns_maps(:,:,:,1,1,1,1,1,1,1,1));
    l1eigs_analytical = sum(Nsamps(:) / N .* squeeze(mdsum(abs(sns_maps).^2,setdiff(1:TIMEMAP_DIM, TIME_DIM))));
end
