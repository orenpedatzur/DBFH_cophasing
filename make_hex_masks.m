function M = make_hex_masks(x_centers, y_centers, hex_radius, image_size)
%MAKE_HEX_MASKS  Build per-spot hexagonal masks stacked in a 3-D array.
%   M = MAKE_HEX_MASKS(x_centers, y_centers, hex_radius, image_size)
%
%   Inputs
%     x_centers, y_centers : vectors of equal length (columns = x, rows = y)
%     hex_radius           : scalar or vector. Distance from center to vertex
%                            (hex circumradius) in pixels.
%     image_size           : [nRows nCols] of the output image(s)
%
%   Output
%     M : logical array of size [nRows nCols N], where N = numel(x_centers).
%         Each "page" M(:,:,k) is a filled, flat-top hexagon centered at
%         (x_centers(k), y_centers(k)) with the specified radius.
%
%   Notes
%     - Orientation is flat-top (horizontal sides). To make pointy-top, set
%       baseAngle = 0 instead of 30 degrees.
%     - Coordinates can be non-integer; poly2mask handles subpixel vertices.
%
%   Example
%     M = make_hex_masks([512], [512], 60, [1024 1024]);
%     imshow(M(:,:,1));

    arguments
        x_centers (:,1) double
        y_centers (:,1) double
        hex_radius (1,:) double {mustBePositive}
        image_size (1,2) double {mustBeInteger, mustBePositive}
    end

    N = numel(x_centers);
    if numel(y_centers) ~= N
        error('x_centers and y_centers must have the same length.');
    end

    % Allow radius per-hex or a single scalar.
    if isscalar(hex_radius)
        hex_radius = repmat(hex_radius, N, 1);
    elseif numel(hex_radius) ~= N
        error('hex_radius must be scalar or length N.');
    else
        hex_radius = hex_radius(:);
    end

    nRows = image_size(1);
    nCols = image_size(2);
    M = false(nRows, nCols, N, 'like', true);

    % Flat-top hexagon: vertices at 30° + k*60°
    baseAngle = 0*deg2rad(30); % 0 for flat top, 30 for pointy top
    ang = baseAngle + (0:5) * (pi/3);   % 6 vertices

    for k = 1:N
        R  = hex_radius(k);
        xc = x_centers(k);
        yc = y_centers(k);

        xv = xc + R * cos(ang);
        yv = yc + R * sin(ang);

        % Rasterize polygon into a binary mask
        if ~isnan(xc*yc)
            mask = poly2mask(xv, yv, nRows, nCols);
            M(:,:,k) = mask;
        else
            M(:,:,k) = zeros(image_size);
        end
    end
end
