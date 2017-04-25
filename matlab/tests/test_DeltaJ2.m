% Check that two methods for computing delta J are equivalent
clear, close all;

ny = 24;
nz = 12;
ncoils = 8;
nmaps = 2;
nt = 7;
ntmaps = 4;

dims = [1 ny nz ncoils nmaps nt ntmaps];

sns_maps = randn(dims) + sqrt(-1)*randn(dims);

w = buildW(sns_maps);

mask = randn(ny,nz,nt) < 0.4;


dJ1 = getDeltaJ(mask, w, 1);
dJ2 = getDeltaJ(mask, w, 0);

threshTest(RelativeError(dJ1, dJ2), -100, 'C++ vs Matlab calculation of Delta J');



