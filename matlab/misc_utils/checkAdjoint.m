% Check that the operator A's adjoint is correctly implemented
% by checking the property
% < Z, A*W> = <A'*Z, W>
% here, Z,W ~ CN(0, sigma^2 I)
% checkAdjoint(A, At, szX, tol)
% INPUTS:
%   A,At    forward operator and adjoint, function handle or operator
%   szX     domain of A
%   tol     Tolerance (e.g. 1e-7)
function checkAdjoint(A, At, szX, tol)
    if nargin < 3
        tol = 1e-6;
    end

    % create data
    W = crandnu(szX);
    if isa(A, 'function_handle')
        codomain = (A(W) ~= 0);
    else
        codomain = (A*W ~= 0);
    end
    szY = size(codomain);
    Z = codomain.*crandnu(szY);
    
    if isa(A, 'function_handle')
        Y = A(W);
    else
        Y = A*W;
    end
    if isa(At, 'function_handle')
        X = At(Z);
    else
        X = At*Z;
    end

    % compare results
    err = abs( vec(Z)'*vec(Y) - vec(X)'*vec(W) );
    threshTest(err, tol, 'Adjoint property');
end

% n = crandnu(sz)
% complex Gaussian noise
function n = crandnu(varargin)
    n = randn(varargin{:}) + sqrt(-1)*randn(varargin{:});
    n = n / sqrt(2);
end
