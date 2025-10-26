function segments = make_segments(img_res, seg_flat_diam_px, focal_length_m)
%MAKE_SEGMENTS  Per-segment coefficients for a diffraction-limited focus.
% segments = make_segments(img_res, seg_flat_diam_px, focal_length)
%
% Outputs a struct 'segments' with fields:
%   pistons : [37x1] piston means (rad)
%   tilt_x  : [37x1] tilt-x means (rad)
%   tilt_y  : [37x1] tilt-y means (rad)
%   defocus : [37x1] defocus means (rad), same value for all segments
%
% Notes
% - Uses the same hex basis as hex_wavefront_random: piston, u (tilt_x),
%   v (tilt_y), and (u^2+v^2) for defocus, where u=x/R, v=y/R and
%   R = seg_flat_diam_px/sqrt(3) (center->vertex) in *pixels*.
% - The defocus coefficient is set from the requested focal_length (meters).
%   You can change lambda/pixel_pitch below to your simulation values.

    %# ---- simulation constants (edit to your setup) ----
    lambda      = 550e-9;     % wavelength [m]
    seg_flat_diam_m = 1; % mirror size [m]
    pixel_pitch = seg_flat_diam_m/seg_flat_diam_px;      % pixel pitch [m/pixel]

    %# ---- hex geometry: flat-top layout ----
    Rpx  = seg_flat_diam_px / sqrt(3);     % [pixels] center->vertex
    aX   = 1.5 * Rpx;                       % axial->pixel (x)
    aY   = sqrt(3) * Rpx;                   % axial->pixel (y)
    axial = generate_axial_37();            % [q r] ring walk used elsewhere

    q = axial(:,1);  r = axial(:,2);
    Xc = aX * q;                            % pixel centers (origin at mosaic center)
    Yc = aY * (r + q/2);
    Uc = Xc / Rpx;                          % normalized centers u0 = x/R
    Vc = Yc / Rpx;                          % normalized centers v0 = y/R

    %# ---- physical defocus -> basis coefficient (radians) ----
    % Global phase: phi(u,v) = c_def*(u^2+v^2), with converging curvature negative.
    Rm   = Rpx * pixel_pitch;               % [m]
    if isinf(focal_length_m) || focal_length_m==0
        c_def = 0;
    else
        c_def = -(pi/lambda) * (Rm^2) * (1/focal_length_m);   % [rad]
    end

    %# ---- set per-segment means for a smooth global paraboloid ----
    defocus = repmat(c_def, 37, 1);
    tilt_x  = 2*c_def * Uc;                 % ∂phi/∂u at (u0,v0)
    tilt_y  = 2*c_def * Vc;                 % ∂phi/∂v at (u0,v0)
    pistons = c_def * (Uc.^2 + Vc.^2);      % phi(u0,v0)

    %# ---- package result ----
    segments = struct( ...
        'pistons', pistons(:), ...
        'tilt_x',  tilt_x(:),  ...
        'tilt_y',  tilt_y(:),  ...
        'defocus', defocus(:), ...
        'meta', struct('img_res', img_res, ...
                       'seg_flat_diam_px', seg_flat_diam_px, ...
                       'focal_length_m', focal_length_m, ...
                       'lambda_m', lambda, ...
                       'pixel_pitch_m', pixel_pitch) );
end
