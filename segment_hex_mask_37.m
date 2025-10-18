function mask = segment_hex_mask_37(indices, img_res, seg_flat_diam_px)
%SEGMENT_HEX_MASK_37  Logical mask for one or more of 37 hex segments (flat-top).
%
%   mask = segment_hex_mask_37(indices, img_res, seg_flat_diam_px)
%     indices: vector of integers in [1..37]  (1=center)
%     img_res: square image size (e.g., 2048)
%     seg_flat_diam_px: flat-to-flat diameter (pixels) of one hex
%
% Geometry (flat-top): flat-to-flat = √3 * R, where R is center→vertex radius.

    arguments
        indices (1,:) {mustBeInteger,mustBePositive,mustBeLessThanOrEqual(indices,37)}
        img_res (1,1) {mustBeInteger,mustBePositive}
        seg_flat_diam_px (1,1) {mustBePositive}
    end

    R  = seg_flat_diam_px / sqrt(3);  % center→vertex radius
    aX = 1.5 * R;                      % axial→pixel (flat-top)
    aY = sqrt(3) * R;

    axial = generate_axial_37();   % [q r] for 37 segments
    cx = (img_res + 1) / 2;
    cy = (img_res + 1) / 2;

    [X, Y] = meshgrid(1:img_res, 1:img_res);
    mask = false(img_res, img_res);

    ang = deg2rad(0 + 60*(0:5));  % flat-top vertices

    for idx = indices
        q = axial(idx,1);  r = axial(idx,2);
        xc = cx + aX*q;
        yc = cy + aY*(r + q/2);

        vx = xc + R*cos(ang);
        vy = yc + R*sin(ang);

        mask = mask | inpolygon(X, Y, vx, vy);
    end
end
