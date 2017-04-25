% Cost increase due to adding a sample
%
% c = getDeltaJ(mask, w, s)
% INPUTS:
%   mask = sampling mask [ny nz]
%   w    = weighting function
%   s    = new sample subscript (3-vector)
% OUTPUTS:
%   c    = cost change from adding s (scalar)
%   dp   = change in p(\Delta k)
function [deltaJ, dp] = getDeltaJ1(mask, w, snew)

% Check inputs
assert(all(mask(:)==0 | mask(:) == 1));

% Mask dimensions
[ny nz nm] = size(mask);

% Compute deltaJ sample add cost
fprintf('Computing add cost\n');

[sy, sz, sm] = ind2sub(size(mask), find(mask(:)~=0));
Ns = sum(mask(:));
ky_new = snew(1);
kz_new = snew(2);
m_new = snew(3);
dp = zeros(ny, nz, nm, nm);

for i = 1:Ns
    dp(mod(sy(i) - ky_new, ny)+1, mod(sz(i) - kz_new, nz)+1,sm(i),m_new,1) = ...
        dp(mod(sy(i) - ky_new, ny)+1, mod(sz(i) - kz_new, nz)+1,sm(i),m_new,1) + 1;
    dp(mod(ky_new - sy(i), ny)+1, mod(kz_new - sz(i), nz)+1,m_new,sm(i),1) = ...
        dp(mod(ky_new - sy(i), ny)+1, mod(kz_new - sz(i), nz)+1,m_new,sm(i),1) + 1;
end
dp(1,1,m_new) = dp(1,1,m_new) + 1;

deltaJ = sum(vec(dp.*w));


end % end getDeltaJ1
