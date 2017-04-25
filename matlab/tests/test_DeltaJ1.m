clear, close all;

if 0
    cd('../src');
    fprintf('Building getDeltaJ_mex.c\n');
    mex -DMEX_COMPILE_FLAG -O CFLAGS="\$CFLAGS -std=c99 -pedantic" getDeltaJ_mex.c dda_utils.c debug.c misc.c multind.c cfl.c
    cd('../tests');
end

ncoils = 3; %3;
nmaps = 2;
nechoes = 3;
ny = 10;
nz = 12;
sns_dims = [1 ny nz ncoils nmaps nechoes];

sns_maps = (randn(sns_dims)  + randn(sns_dims)*sqrt(-1))/sqrt(2);

w = buildW(sns_maps);

mask = rand([ny nz nechoes]) < 0.2;

deltaJ = getDeltaJ(mask, w, 0);

p = get_dd(mask);

cost = sum(w(:).*conj(p(:)));

[sy sz sm] = ind2sub(size(mask), find(mask(:)==0));

sprime = ceil(rand*length(sy));
sy1 = sy(sprime);
sz1 = sz(sprime);
sm1 = sm(sprime);

%[deltaJ1, dp] = getDeltaJ1(mask, w, [sy1, sz1, sm1]);

mask2 = mask;
mask2(sy1,sz1,sm1) = 1;

p2 = get_dd(mask2);

cost2 = sum(w(:).*conj(p2(:)));

cost2_est = cost + deltaJ(sy1,sz1,sm1);

threshTest(RelativeError(cost2, cost2_est), -100, 'C + Delta C');





