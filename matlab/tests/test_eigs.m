% Test the conditioning metrics on 1D sense problems
% Verify relationships in the paper for general 
% (multi-map, multi-coil, multi-echo) block matrices
clear, close all;
% --- Settings
MAPS_DIM = 5;
MAPS2_DIM = 7;
TE_DIM = 6;
COIL_DIM = 4;
MAPS2_DIM = TE_DIM+1;
nMonte = 1000; %10000;
N = 10;
sns_type = 0; % 0 => random
random_sampling = true;
R = 2.3;
nt = 4;
ncoils = 6;
nmaps = 5;
ntmaps = 3;
% ----

% Fourier
Fmtx = kron(eye(nt*ncoils), 1/sqrt(N)*dftmtx(N));  

img_dims = ones(1,8);
img_dims([1:3, MAPS2_DIM, MAPS_DIM]) = [1 1 N ntmaps nmaps];
Sx_dims = ones(1,8);
Sx_dims([1:3, TE_DIM, COIL_DIM]) = [1 1 N nt ncoils];

% Sensitivities
if sns_type == 0 
    disp('Random sensitivities');
    sns_maps = randn(1,1,N,ncoils,nmaps,nt,ntmaps) + randn(1,1,N,ncoils,nmaps,nt,ntmaps)*sqrt(-1);
    % S operator forward
    Sop = @(X) mdsum(sns_maps.*mdrepmat(X, [COIL_DIM TE_DIM], [ncoils nt]), [MAPS_DIM, MAPS2_DIM]);
    % S operator adjoint
    Sopt = @(X) mdsum(conj(sns_maps).*mdrepmat(X, [MAPS_DIM, MAPS2_DIM], [nmaps, ntmaps]), [COIL_DIM, TE_DIM]);

    checkAdjoint(Sop, Sopt, img_dims, 1e-7);

    Smtx = zeros(prod(Sx_dims), prod(img_dims));
    for i = 1:prod(img_dims)
        %ei = zeros(prod(img_dims),1);
        ei = zeros(img_dims);
        ei(i) = 1;
        Smtx(:,i) = vec(Sop(ei));
    end
end

Mmtx = [];
for t = 1:nt
    m_t = rand(N,1) < 1/R;
    mask(:,:,t) = m_t;
for c = 1:ncoils
    Mmtx(end+1:end+N,end+1:end+N) = diag(m_t);
end
end
for t = 1:nt
    Reff(t) = numel(mask(:,:,t)) / sum(vec(mask(:,:,t)));
    Nt(t) = sum(vec(mask(:,:,t)));
end

% Key quantities
E = Mmtx * Fmtx * Smtx;
wmtx = buildW(sns_maps);
p = get_dd(mask);

l1eigs_analytical = get_l1eigs(sns_maps, Nt); 
l2eigs_analytical = sum(wmtx(:).*conj(p(:)));
%l2eigs_analytical  / (1/N^2 * l1eigs_analytical^2);

% Compute some statistics
disp('Eigendecomposition');
eigs_EhE = eig(E'*E);
l1eigs_EhE      = norm(eigs_EhE,1);
l2eigs_EhE      = norm(eigs_EhE,2)^2;

% Estimate norms by simulation
fprintf('Estimatem norms by simulation...');
An1 = NaN*ones(nMonte,1);
An2 = NaN*ones(nMonte,1);
for k = 1:nMonte
    noise = 1/sqrt(2) * (randn(size(E,2),1) + sqrt(-1)*randn(size(E,2),1));
    An1(k) = norm(E*noise,2)^2;
    An2(k) = norm(E'*E*noise,2)^2;
end
fprintf('Done!\n');

% Print results to check formulas
fprintf('---- Sum of eigenvalues  ---- \n');
fprintf('E|En|^2:\t%f +/- %.2f\n', mean(An1), std(An1)/sqrt(numel(An1)));
fprintf('|E|F^2:\t%f\n', norm(E, 'fro')^2);
fprintf('sum lambda_k(E''E):\t%f\n', l1eigs_EhE);
fprintf('Analytical:\t%f\n', l1eigs_analytical);

fprintf('---- Sum of squares of the eigenvalues  ---- \n');
fprintf('sum lambda_k(E''E)^2:\t%f\n', l2eigs_EhE);
fprintf('E|E''E*E''E*n|^2:    \t%f +/- %.2f\n', mean(An2), std(An2)/sqrt(numel(An2)));
fprintf('|E''E|F^2:           \t%f\n', norm(E'*E, 'fro')^2);
fprintf('<w, p>:              \t%f\n', l2eigs_analytical);


