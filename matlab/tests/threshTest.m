% threshTest(err, thresh, [testname])
%   print pass or fail for a test where we must have err < thresh
function threshTest(err, thresh, testname)
    if nargin == 3
        fprintf('---- %s ----\n', testname);
    end

    if err < thresh
        fprintf('PASSED\terr: %f < %f\n', err, thresh);
    else
        fprintf('FAILED!\terr: %f > %f\n', err, thresh);
    end
end

