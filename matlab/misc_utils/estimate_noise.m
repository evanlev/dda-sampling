% estimate_noise(mask, AtAmtx, w, nMonte, name)
%
function [eps2_mean] = estimate_noise(AtAmtx, w, nMonte, name)

    img_dims = size(mask);
    N = prod(img_dims);

    % Empirical method
    for k = 1:nMonte
        noise = randn(img_dims);
        eps2(k) = norm(vec(AtAmtx*noise))^2 / N; % or ||A noise||^2/ ||noise||^2
    end
    eps2_mean = mean(eps2);
    conf = 1.96*std(eps2)/sqrt(nMonte);
    
    % Print result
    fprintf('%s:\t Noise level for EtE*eps: %f +/- %f(empirical)\n', ...
            name, mean(eps2), conf);
end
