clear, close all;
#addpath('../utils');
TE_DIM = 6;
TEMAP_DIM = 7;
COIL_DIM = 4;
disp('Approximate best candidate sampling demo...');
% ---- Settings ------
% Sampling optimizer
show_ebc = false; % set to true to show exact best candidate sampling
show_gfactor = false; % Set true to show g-factor
tol = 0;
R = 4;
K = 32;               % Support of w to use
nMonte = 200; 
nMonteFull  = 400; 
bcfast_path = '../../../cpp/bcfast/bcfast'
% -------------------

%% Set up save file
savestr = sprintf('results/res_approx_R%d.mat', R);
fprintf('Data will be saved to %s\n', savestr);

%% Set up experiments
% Coil sensitivities
load('../../data/vfa4d_ccomp_sl130.mat');

if ~exist(bcfast_path, 'file')
    error(["Please compile the bcfast executable and set the path to it. Program does not exist at " bcfast_path])
end

if isempty(getenv('TOOLBOX_PATH'))
		error('This demo requires BART to be installed and the environment variable TOOLBOX_PATH to be set.');
end

% Synthetic data
f_coil_combine = @(X) bart_0700('pics -w1 -S -d0 -i1 ', X, sns_maps); 
f_sense =        @(X) bart_0700('pics -w1 -S -d0 -i75', X, sns_maps);

%% Noise map full
disp('g-factor map for fully sampled image');
f_rms = @(X) sqrt(mean(X(:).^2));
sigma = f_rms(kData)*1e-1;
kData = (randn(size(kData))*sqrt(-1) + randn(size(kData))) * sigma * 100/sqrt(2);

if show_gfactor
    [gMapFull, reconFull] = gfactor(ones(size(kData)), nMonteFull, 1, tol, f_coil_combine, 1);
    supp = gMapFull ~= 0;
    gMapFull(~supp) = Inf;
end

%% Create masks
masks = [];
names = {};
% -- Poisson-disc
mask_pd = bart_0700(sprintf('poisson -Y%d -Z%d -y%f -z%f', size(kData,2), size(kData,3), ...
        sqrt(R), sqrt(R)));
masks = cat(4,masks,mask_pd);
names{end+1} = 'Poisson-disc';
maxSamples = sum(mask_pd(:));

%% Get w 
wmtx = buildW(sns_maps);

% -- Approximate best candidate
tic;
mask_abc = best_candidate(bcfast_path, wmtx, maxSamples, 'K', K);
fprintf('Approximate best candidate with K = %d took %f sec\n', K, toc);
mask_abc = shiftdim(mask_abc, -1);
masks = cat(4,masks,mask_abc);
names{end+1} = 'Approximate BC';

% -- Best candidate 
if show_ebc
    tic;
    mask_ebc = best_candidate(bcfast_path, wmtx, maxSamples);
    t0 = toc;
    fprintf('Exact best candidate sampling took %f sec\n', t0);
    mask_ebc = shiftdim(mask_ebc, -1);
    masks = cat(4,masks,mask_ebc);
    names{end+1} = 'Exact BC';
end

%% Differential distributions for all masks
fprintf('Get differential distributions...\n');
for k = 1:size(masks,4)
    diff_dist(:,:,:,k) = shiftdim(get_dd(shiftdim(masks(1,:,:,k),1)),-1);
end

%% <w,p> cost
fprintf('Pattern     \t<w,p>\n');
l1eigs_analytical = 1/R*sum(abs(sns_maps(:)).^2);
wmtx = wmtx / (1/numel(sns_maps(:,:,:,1))*l1eigs_analytical^2);
for k = 1:size(diff_dist,4)
    wp_cost(k) = sum(wmtx(:).*vec(diff_dist(:,:,:,k)));
    fprintf('%s\t%.1f\n', names{k}, wp_cost(k));
end

%% G-factor maps
if show_gfactor
    fprintf('g-factor: WARNING slow...\n');
    for k = 1:size(masks,4)
        fprintf('g-factor (%d/%d)\n', k, size(masks,4));
        % g-factor maps
        [gMap(:,:,:,k), im{k}] = gfactor(mdrepmat(masks(:,:,:,k), COIL_DIM, size(sns_maps,COIL_DIM)), nMonte, 1, tol, f_sense, 1);
        
        % g-factor stats
        gMap(:,:,:,k) = gMap(:,:,:,k) ./ gMapFull;
        
        gRMS(k) = f_rms(vec(gMap(:,:,:,k)));
        gMax(k) = max(vec(gMap(:,:,:,k)));
        gMean(k) = mean(vec(gMap(:,:,:,k)));
    end
end

%% Cost in k-space
fprintf('Get Delta J...\n');
for k = 1:size(masks,4)
    deltaJ(:,:,k) = getDeltaJ(squeeze(masks(:,:,:,k)), wmtx);
end
fprintf('Saving %s...', savestr);

if ~exist('OCTAVE_VERSION', 'builtin') ~= 0
    save(savestr, '-v7.3');
end

fprintf('%s done\n', mfilename('fullpath'));

