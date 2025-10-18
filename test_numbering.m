function test_tumbering(img_res,seg_flat_diam_px)

figure; clf;
hold on;

for i = 1:37
    % Get mask for segment i
    M = segment_hex_mask_37(i, img_res, seg_flat_diam_px);

    % Find its center
    [y, x] = find(M);
    xc = mean(x);
    yc = mean(y);

    % Show mask outline
    B = bwboundaries(M);
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1),'color', 0.5*[1 1 1], 'LineWidth', 2);
    end

    % Label with segment number
    text(xc+cosd(60)*seg_flat_diam_px/2.4, yc-sind(60)*seg_flat_diam_px/2.4, num2str(i), ...
        'FontWeight', 'bold', 'FontSize', 16, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end

axis equal ij;
xlim([1 img_res]);
ylim([1 img_res]);
axis off;
title('Segment IDs for 37-segment primary mirror');
end