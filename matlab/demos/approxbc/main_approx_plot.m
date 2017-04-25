% Plot
clear, close all;
load('results/res_approx_R4.mat');
show_sk = true;
gamma_p = 0.2;
gamma_psf = 0.2;
h = {}; hc = {};

% Permute arrays
if show_gfactor
    gMap = permute(gMap, [1 3 2 4]);
end
masks = permute(masks, [1 3 2 4 5]);
wmtx = permute(wmtx, [2 1]);
diff_dist = permute(diff_dist, [1 3 2 4]);
sns_maps = permute(sns_maps, [1 3 2 4]);
deltaJ = abs(permute(deltaJ, [2 1 3]));
npats = size(masks,4);
clear gMapFull;
clear kData;

% w(Delta k)
N = 30;
figure;
w_disp = mdcrop(fftshift(fftshift(wmtx, 1), 2), 1:2, [N N]);
surf(w_disp);
title('w(\Delta k)');

% PSFs
count = 1;
psfs = mdfft(masks, [2 3], 1, 1, false);
figure;
for k = 1:size(masks,4)
    subplot(3,npats,count);
    count = count + 1;
    psf_disp = fftshift(fftshift(shiftdim(abs(psfs(:,:,:,k)),1), 1), 2).^gamma_psf;
    psf_disp(end/2+1,end/2+1) = psf_disp(end/2+1,end/2+1) / 4;
    %psf_disp(end/2 + (-1:1),end/2 + (-1:1)) = psf_disp(end/2 + (-1:1),end/2 + (-1:1))/6;
    imagesc(psf_disp, [0, max(psf_disp(:))]);
    colormap gray;
    title(names{k});
end


surf2 = @(X) surf(X, 'FaceAlpha', 1);
% p( Delta k)
N = 30;
for k = 1:size(diff_dist,4)
    subplot(3,npats,count);
    count = count + 1;
    p_disp = mdcrop(fftshift(fftshift(shiftdim(diff_dist(:,:,:,k), 1), 1), 2), 1:2, [N N]).^gamma_p;
    imagesc(p_disp, [0, max(p_disp(:))]);
    colormap gray;
    title(['p(\Delta k): ', names{k}]);
end

% s(k)
N = 14;
deltaJ_seg = mdcrop(deltaJ, 1:2, [N N]);
maxDeltaJ = max(deltaJ_seg(:));
for k = 1:size(masks,4)
    subplot(3,npats,count);
    count = count + 1;
    imagesc(deltaJ_seg(:,:,k), [0, maxDeltaJ]);
    colormap gray;
    hold on;
    mask_seg = shiftdim(mdcrop(masks(1,:,:,k), [2 3], [N N]),1);
    [sz, sy] = ind2sub(size(mask_seg), find(mask_seg(:)));
    plot(sy(:)-1, sz(:)-1, 'k.', 'MarkerSize', 1);
    hc{end+1} = colorbar;
    title(names{k});
end

% g-factor
if show_gfactor
    yl = [0, max(vec(gMap(:,:,:,1)))];
    figure;
    for k = 1:npats
        subplot(1,npats,k);
        imagesc(gMap(:,:,:,k), yl);
        title(names{k});
        colormap jet;
        hc{end+1} = colorbar;
    end
end


