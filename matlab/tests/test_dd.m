clear, close all;

ny = 32;
nz = 16;
nt = 4;

mask = rand(ny,nz,nt) < 4;

p1 = get_dd(mask, 'Complex', 1);
p2 = get_dd(mask, 'Complex', 0); 

threshTest(RelativeError(p1,p2), -80, 'BART vs MATLAB p(Delta k) calculation');






