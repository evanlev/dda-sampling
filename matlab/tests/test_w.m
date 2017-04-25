clear, close all;

ny = 24;
nz = 12;
ncoils = 8;
nmaps = 2;
nt = 7;
ntmaps = 4;

dims = [1 ny nz ncoils nmaps nt ntmaps];

sns_maps = randn(dims) + sqrt(-1)*randn(dims);

w1 = buildW(sns_maps, true);
w2 = buildW(sns_maps, false);
fprintf('\n');

threshTest(RelativeError(w1, w2), -100, 'BART vs MATLAB calculation of w');



