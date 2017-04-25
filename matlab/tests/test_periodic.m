clear, close all;

% Test the periodicity property of p
N0 = 10;
N = N0*10;

% Generate random periodic pattern
s0 = rand(N0,N0) < 0.4;
s = repmat(s0, [N/N0 N/N0]);

p2 = get_dd(s);
p1 = repmat(get_dd(s0), [N/N0, N/N0])*(N/N0)^2;

threshTest(RelativeError(p1,p2), -80, 'Periodicity of p');



