clear, close all;
R = 6;

% Generate CAIPIRINHA patterns and try to predict the noise
% Compare this to the RMS g-factor

% ---- Settings ------
caip_sampling = true;   % False => all periodic patterns
nMonte = 2;
% --------------------
fprintf('R = %d, caip_sampling = %d\n', R, caip_sampling);

% --- Set up data
load('../../data/vfa4d_ccomp_sl130.mat');
[nx ny nz ncoils] = size(kData);
% Pad so that the dimensions are a multiple of the period
ny = ny + R - mod(ny, R);
nz = nz + R - mod(nz, R);
N = ny*nz;
kData = zpad(kData, [nx, ny nz ncoils]);
sns_maps = bart_0209('ecalib -m1 -c0.4', kData);

% Make weight vector: <W, p(deltak)> = || E'E epsilon ||_2^2, epsilon ~ N(0,I)
w = buildW(sns_maps);

%% --- Test CAIPIRINHA masks
if caip_sampling
    % Only CAIP kernels
    cells = getKernels(R, 'CAIPIRINHA');
else
    % Arbitrary kernels
    cells = getKernels(R, 'All');
    % Shuffle 
    cells = cells(:,:,randperm(size(cells,3)));
    % Limit to 50 
    cells = cells(:,:,1:min(50, size(cells,3)));
end
npats = size(cells,3);

% Reconstruction functions 
f_coil_combine = @(X) bart_0209('pics -S -d0 -i1 -w1', X, sns_maps);
f_sense = @(X) bart_0209('pics -S -d0 -w1', X, sns_maps);

% Noise to add 
sigma = norm(kData(:))/numel(kData(:))*9e-2;

% statistics and g-factor
gMap = zeros(nx,ny,nz,npats);

fprintf('Noise estimation for periodic sampling\n');
for k = 1:npats
    fprintf('Kernel %d/%d\n', k, npats);
    cell = cells(:,:,k); 
    repgrid = ones(ceil(ny/size(cell,1)), ceil(nz/size(cell,2)));
    mask_ca = kron(repgrid, cell);

    masks(:,:,:,k) = shiftdim(mask_ca, -1);

    diff_dist(:,:,k) = kron(repgrid, get_dd(cell)*N/(R^2));

    wp(k) = sum(vec(w.*diff_dist(:,:,k)));
    wp_normalized(k) = wp(k) / (1/N * get_l1eigs(sns_maps, sum(vec(masks(:,:,:,1))))^2);

    gk = gfactor_caip(sns_maps, cell);
    gMap(:,:,:,k) = gk;

    gMean(k) = mean(gk(gk(:)~=0));
    gMax(k) = max(abs(gk(gk(:)~=0)));
    gNorm(k) = sqrt(sum(abs(gk(gk(:)~=0)).^2));
    fprintf('Mean g: %f\n', gMean(k));
    if( gMean(k) == 0 )
        error('g-factor has mean 0, missing -w1 option?');
    end

    deltaJ(:,:,k) = getDeltaJ(mask_ca, w);
end

% --- Save
fname = sprintf('results/res_caip_R%d', R);
if( ~caip_sampling )
    fname = [fname, '_arb'];
end
save(fname, '-v7.3');

main_caip_plot;

% --------------------------
fprintf('Done main_caip.m\n');





