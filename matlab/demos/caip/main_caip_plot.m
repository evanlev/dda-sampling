%% ---------- Make some plots ----------------------
close all;
h = {};
hc = {};

high_g = gMax > 5*median(gMax(:));
wp_normalized2 = wp_normalized(~high_g);
gMap2 = gMap(:,:,:,~high_g);
deltaJ2 = deltaJ(:,:,~high_g);
npats = size(gMap2, 4);


% --- Figure 1: sampling patterns and delta J
N = 1.5*R;
deltaJ_seg = mdcrop(deltaJ2, 1:2, [N N]);
mindeltaJ = min(deltaJ_seg(:));
maxdeltaJ = max(deltaJ_seg(:));
h{end+1} = figure;
for k = 1:npats
    %h = plot_sk(masks(:,:,:,k), N, deltaJ
    if npats <= 6
        subplot(2,3,k);
    elseif npats <= 9
        subplot(3,3,k);
    elseif npats <= 12
        subplot(3,4,k);
    else
        subplot(5,5,k);
    end
    imagesc(deltaJ_seg(:,:,k), [mindeltaJ, maxdeltaJ]);
    colormap gray;
    hold on;
    mask_seg = shiftdim(mdcrop(masks(1,:,:,k), [2 3], [N N]),1);
    [sz, sy] = ind2sub(size(mask_seg), find(mask_seg(:)));
    plot(sy(:)-1, sz(:)-1, 'k.');
    title(sprintf('<w,p>: %.2f', wp_normalized(k)));
    %hc{end+1} = colorbar;
end

% --- g-factor
h{end+1} = figure;

yl = [0 max(gMap2(:))];
for k = 1:npats
    if npats <= 6
        subplot(2,3,k);
    elseif npats <= 9
        subplot(3,3,k);
    elseif npats <= 12
        subplot(3,4,k);
    else
        subplot(5,5,k);
    end
    imagesc(shiftdim(gMap2(:,:,:,k),1), yl);
    axis image
    axis tight;
    title(sprintf('<w,p>: %.2f', wp_normalized(k)));
    colormap jet;
    hc{end+1} = colorbar;
end

