function unstacked_segments = unstack_segment_tilts(segments,unstack_ratio)
%UNSTACK_SEGMENT_TILTS  Remove nominal focusing tilts from each segment.
%   unstacked = unstack_segment_tilts(segments)
%
% Subtracts the per-segment tilt means that come from the global quadratic
% phase matching the focal length in segments.meta.focal_length_m:
%   phi(u,v) = c0*(u^2+v^2),  c0 = -(pi/lambda)*(R_m^2/f0)
% Nominal tilts are grad(phi) at each segment center (u0,v0):
%   tilt_x_nom = 2*c0*u0,   tilt_y_nom = 2*c0*v0
%
% Only tilt_x and tilt_y are adjusted. Pistons/defocus/meta are unchanged.

    % ---- checks ----
    need = {'img_res','seg_flat_diam_px','focal_length_m','lambda_m','pixel_pitch_m'};
    assert(isfield(segments,'meta') && all(isfield(segments.meta, need)), ...
        'segments.meta must contain: %s', strjoin(need, ', '));

    N = numel(segments.tilt_x);
    assert(N==37, 'segments fields must have 37 elements.');

    % preserve original shape (row/col)
    shp = size(segments.tilt_x);

    % ---- hex geometry & segment centers in normalized coords (u0,v0) ----
    Rpx  = segments.meta.seg_flat_diam_px / sqrt(3);  % center->vertex [px]
    aX   = 1.5 * Rpx;
    aY   = sqrt(3) * Rpx;
    % aY   = 1.5 * Rpx;
    % aX   = sqrt(3) * Rpx;
    axial = generate_axial_37();      % [q r] layout consistent with prior code
    q = axial(:,1);  r = axial(:,2);

    Xc = aX * q;                      % centers in pixels (mosaic-centered)
    Yc = aY * (r + q/2);

    u0 = Xc / Rpx;                    % normalized centers
    v0 = Yc / Rpx;

    % ---- nominal focusing coefficient from focal length ----
    f0   = segments.meta.focal_length_m;
    lam  = segments.meta.lambda_m;
    pixm = segments.meta.pixel_pitch_m;

    R_m  = Rpx * pixm;                                   % [m]
    if isinf(f0) || f0==0
        c0 = 0;
    else
        c0 = -(pi/lam) * (R_m^2) * (1/f0);               % [rad]
    end

    tx_nom = 2*c0*u0;
    ty_nom = 2*c0*v0;

    % ---- build output (subtract nominal tilts) ----
    unstacked_segments       = segments;         % copy everything
    unstacked_segments.tilt_x = reshape(segments.tilt_x(:) - unstack_ratio*tx_nom, shp);
    unstacked_segments.tilt_y = reshape(segments.tilt_y(:) - unstack_ratio*ty_nom, shp);
    % pistons/defocus/meta unchanged
end
