function [U, phase, mask, coeffs, terms] = hex_wavefront_random( ...
    seg_flat_diam_px, terms, mu_sigma, seed, center_offset)
%HEX_WAVEFRONT_RANDOM  Random-coefficient wavefront over a flat-top hex tile.
%
%   [U, phase, mask, coeffs, terms] = hex_wavefront_random(seg_flat_diam_px, terms, mu_sigma, seed)
%   [U, phase, mask, coeffs, terms] = hex_wavefront_random(..., center_offset)
%
% Inputs
%   seg_flat_diam_px : flat-to-flat size of the hex in *pixels*.
%   terms            : cellstr of basis terms. Supported:
%                      'piston','tilt_x','tilt_y','defocus',
%                      'x2','xy','y2','x3','x2y','xy2','y3','x4','x3y','x2y2','xy3','y4'
%   mu_sigma         : [N x 2] array of [mean, std] (radians) per term in 'terms'.
%   seed             : (optional) integer RNG seed for reproducibility.
%   center_offset    : (optional) [dx,dy] sub-pixel offset in *pixels* to shift the
%                      tile's center relative to its array center. Default [0 0].
%
% Outputs
%   U      : complex wavefront, U = mask .* exp(1i * phase)
%   phase  : phase map (radians), same size as mask
%   mask   : logical hex mask (flat-top), size ~ [2R x seg_flat_diam_px] (odd×odd)
%   coeffs : drawn coefficients (1 x N) in radians
%   terms  : echo of the requested terms (cellstr)
%
% Notes
% - Normalized coords: u = x/R, v = y/R, with R = seg_flat_diam_px/sqrt(3) (center→vertex).
% - For perfect registration with external mosaics, pass the fractional center
%   via center_offset. If omitted, the tile is centered on its own pixel grid.

    arguments
        seg_flat_diam_px (1,1) {mustBePositive}
        terms cell {mustBeVector}
        mu_sigma double {mustBeReal, mustBeFinite}
        seed double {mustBeScalarOrEmpty} = []
        center_offset (1,2) double {mustBeReal} = [0 0]
    end

    if size(mu_sigma,1) ~= numel(terms) || size(mu_sigma,2) ~= 2
        error('mu_sigma must be [numel(terms) x 2] of [mean std].');
    end

    if ~isempty(seed)
        rng(seed);
    end

    % ---- Hex geometry (flat-top) and pixel canvas ----
    % R = center-to-vertex radius (pixels); flat-to-flat = sqrt(3)*R
    R  = seg_flat_diam_px / sqrt(3);

    % Bounding box: width = 2R (x), height = seg_flat_diam_px (y)
    nx = round(2*R);
    ny = round(seg_flat_diam_px);

    % Force odd sizes so the array center is exactly a pixel
    if mod(nx,2)==0, nx = nx + 1; end
    if mod(ny,2)==0, ny = ny + 1; end

    % Pixel centers with (1,1) at top-left; array center at (cx,cy)
    cx = (nx + 1)/2;
    cy = (ny + 1)/2;

    [X, Y] = meshgrid(1:nx, 1:ny);

    % Apply optional sub-pixel center offset (dx,dy) in pixels
    dx = center_offset(1);
    dy = center_offset(2);
    x = (X - cx) - dx;   % pixels
    y = (Y - cy) - dy;   % pixels

    % ---- Build the flat-top hex mask via polygon ----
    ang = deg2rad(0 + 60*(0:5));      % flat-top vertices
    vx = R*cos(ang);
    vy = R*sin(ang);
    mask = inpolygon(x, y, vx, vy);   % inside or on boundary

    % ---- Normalized coordinates (dimensionless) ----
    u = x / R;    % vertices ~ radius 1
    v = y / R;

    % ---- Assemble phase from requested basis terms ----
    coeffs = zeros(1, numel(terms));
    phase  = zeros(ny, nx);

    for k = 1:numel(terms)
        t  = lower(string(terms{k}));
        mu = mu_sigma(k,1);
        sg = mu_sigma(k,2);
        c  = mu + sg*randn(1);       % coefficient in radians
        coeffs(k) = c;

        B = basis_term(t, u, v);     % basis map for this term
        phase = phase + c * B;
    end

    % Zero phase outside the hex (cosmetic)
    phase = phase .* mask;

    % Complex wavefront
    U = mask .* exp(1i * phase);
end

% ---- basis term factory ----
function B = basis_term(name, u, v)
    switch name
        case "piston",  B = ones(size(u));
        case "tilt_x",  B = u;
        case "tilt_y",  B = v;
        case "defocus", B = (u.^2 + v.^2);

        % Quadratics
        case "x2",      B = u.^2;
        case "xy",      B = u.*v;
        case "y2",      B = v.^2;

        % Cubics
        case "x3",      B = u.^3;
        case "x2y",     B = (u.^2).*v;
        case "xy2",     B = u.*(v.^2);
        case "y3",      B = v.^3;

        % Quartics
        case "x4",      B = u.^4;
        case "x3y",     B = (u.^3).*v;
        case "x2y2",    B = (u.^2).*(v.^2);
        case "xy3",     B = u.*(v.^3);
        case "y4",      B = v.^4;

        otherwise
            error('Unknown term "%s". See help for supported names.', name);
    end
end
