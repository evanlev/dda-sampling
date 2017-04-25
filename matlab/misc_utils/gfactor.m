% [gMap, samp_x, mu, sig] = gfactor(ksp_mask, nMonte, sdFull, tol, fptr)
%
% Estimates a g-factor map using a Monte Carlo simulation
% Use noise as input data
%
% INPUTS: 
%   ksp_mask = mask, 1 at sample locations, 0 otherwise. Include coil dimension
%   nMonte   = max # monte carlo repetitions (e.g. 100)
%   sdFull   = compute g-factor as sigma ./ (sqrt(R) * sdFull), 
%              where sigma is the pointwise standard deviation of 
%              all reconstructions
%   tol      = Limit nMonte if standard dev of g-factor estimate is < tol
%   fptr     = recon function, image = fptr(data)
%   sigma    = noise covariance matrix. Must have
%              size(ksp_mask,4) == size(sigma,1), 4 is the coil dim
% OUTPUTS:
%   gMap    = g-factor map
%   samp_x  = sample reconstruction
%   mean    = mean map
%   sigm    = standard deviation map
function [gMap, samp_x, mu, sig] = gfactor(ksp_mask, nMonte, sdFull, tol, fptr, sigma)

if nargin < 6
    sigma = 1; % noise to add
end
% ---- Constants
MIN_RUNS = 20; % min # runs before testing stopping criterion
% --------------

support_mask = sdFull ~= 0; % Avoid dividing by zero

% Start monte carlo loop for recon experiment
assert(all(ksp_mask(:) == 0 | ksp_mask(:) == 1));
        
fprintf('g-factor: %d repetitions', nMonte);

for run = 1:nMonte
    if isequal(sigma, 1)
        noise_data = ksp_mask.*(sigma*crandnu(size(ksp_mask)));
    else
        noise_data = ksp_mask.*mvcrandnu_multicoil(size(ksp_mask), sigma);
    end
    samp_x = fptr(noise_data);
    if ~exist('mu', 'var')
        mu = samp_x / nMonte;
        ss = abs(samp_x).^2 / nMonte;
    else
        mu = mu + samp_x / nMonte;
        ss = ss + abs(samp_x).^2 / nMonte;
    end
    sig = support_mask .* sqrt(ss - abs(mu).^2);
    if 0
        fprintf('g-factor stdev of estimate: %f (min) %f (max) %f (mean)\n', ...
                min(sig(:))/sqrt(run), max(sig(:))/ sqrt(run), ...
                mean(sig(:))/sqrt(run));
    end
    if run > MIN_RUNS && max(sig(:)/sqrt(run)) < tol
        fprintf('breaking after %d / %d\n', run, nMonte);
        break;
    end
    % Uncomment to see all realazations
    %res(:,run) = samp_x(:);
end

% sigma^2 = E[|X|^2] - |E[X]|^2 = M2 - M1.^2
R = numel(ksp_mask) / sum(ksp_mask(:)~=0);
gMap = reshape(sig, size(samp_x)) ./ sqrt(R);

fprintf('g-factor stdev of estimate: %f (min) %f (max) %f (mean)\n', ...
        min(sig(:))/sqrt(run), max(sig(:))/ sqrt(run), ...
        mean(sig(:))/sqrt(run));

end

% n = crandnu(sz), complex white Gaussian noise with unit variance
function n = crandnu(sz)
    n = (sqrt(-1)*randn(sz) + randn(sz))/sqrt(2);
end

% Generate noise for data of size sz with covariance sigma
% Assume sz(4) = # coils
function R = mvcrandnu_multicoil(sz, sigma)
    COIL_DIM = 4;
    assert(ndims(sigma) <= 2 && size(sigma,1) == size(sigma,2));

    mu = zeros(size(sigma,1),1);
    R = 1/sqrt(2) * mvnrnd(mu, sigma, prod(sz)/sz(COIL_DIM)) + ...
        sqrt(-1)/sqrt(2) * mvnrnd(mu, sigma, prod(sz)/sz(COIL_DIM));
    N = size(R,1); 
    D = size(R,2); % # coils
    D1 = prod(sz(1:(COIL_DIM-1)));     % dims 1..3
    D2 = prod(sz((COIL_DIM+1):end));   % dims 5..end
    R = permute(reshape(R, [D1 D2 D]), [1 3 2]); % Noise permuted to size D1 D D2
    R = reshape(R, sz);
end


