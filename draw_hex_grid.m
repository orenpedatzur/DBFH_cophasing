function hex_grid = draw_hex_grid(segments, varargin)
%DRAW_HEX_GRID  Draw hexagon outlines on the simulated pupil canvas.
%   hex_grid = draw_hex_grid(segments)
%   hex_grid = draw_hex_grid(segments,'Thickness',2,'Show',true,'Color',[1 1 0])
%
% Returns:
%   hex_grid  -> HxW logical, true on the hex outlines (dilated perimeters)
%
% Options:
%   'Thickness' : edge thickness in pixels (default 2)
%   'Show'      : if true, display overlay (default true)
%   'Color'     : RGB color for overlay lines (default [1 1 0], yellow)

    % ---- options ----
    p = inputParser;
    p.addParameter('Thickness', 2, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
    p.addParameter('Show', true, @(x)islogical(x)||isscalar(x));
    p.addParameter('Color', [1 1 0], @(c)isnumeric(c)&&numel(c)==3);
    p.parse(varargin{:});
    thick = round(p.Results.Thickness);
    doShow = logical(p.Results.Show);
    col    = p.Results.Color;

    % ---- canvas / geometry ----
    H   = segments.meta.img_res;
    W   = segments.meta.img_res;
    Rpx = segments.meta.seg_flat_diam_px / sqrt(3);     % center->vertex (px)
    aX  = 1.5 * Rpx;                                    % axial spacing x (px)
    aY  = sqrt(3) * Rpx;                                % axial spacing y (px)

    axial = generate_axial_37();                        % [q r] per tile
    q = axial(:,1);  r = axial(:,2);

    % global array center (keep it exactly like your pupil)
    cxG = (W+1)/2;
    cyG = (H+1)/2;

    % center positions (px)
    Xc = cxG + aX*q;
    Yc = cyG + aY*(r + q/2);

    % hexagon vertices (flat-top) in local tile coords
    ang = deg2rad(0 + 60*(0:5));                        % 0,60,...,300 deg
    vx  = Rpx*cos(ang);                                 % local x
    vy  = Rpx*sin(ang);                                 % local y

    % ---- rasterize perimeters ----
    hex_grid = false(H,W);
    for i = 1:numel(Xc)
        % polygon in global pixel coords
        xg = vx + Xc(i);
        yg = vy + Yc(i);

        % filled hex mask
        fill_i = poly2mask(xg, yg, H, W);

        % perimeter pixels
        per_i = bwperim(fill_i);

        % thicken edges (no IPT fallback if imdilate missing)
        try
            per_i = imdilate(per_i, strel('disk', max(1,thick-1)));
        catch
            k = ones(max(1,thick), max(1,thick));
            per_i = conv2(double(per_i), k, 'same') > 0;
        end

        hex_grid = hex_grid | per_i;
    end

    % ---- optional display ----
    if doShow
        % a neutral background
        bg = zeros(H,W);
        imshow(bg,[],'Border','tight'); hold on;
        % draw outlines
        [ry,rx] = find(hex_grid);
        plot(rx, ry, '.', 'Color', col, 'MarkerSize', 1);
        axis image ij; title('Hex grid','Color','w');
        hold off;
    end
end
