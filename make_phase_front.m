function [phi_full, mask_full, U_full, phi_stripped, mask_stripped, U_stripped] = ...
    make_phase_front(segments)
%MAKE_PHASE_FRONT  Render exit-pupil phase from per-segment commands, and
%                  also return the phase with the global parabolic curvature
%                  (for the given focal length) removed.
%
% Inputs
%   segments : struct with fields (each length 37)
%       pistons, tilt_x, tilt_y, defocus    % coefficients in radians
%     and meta sub-struct with:
%       meta.img_res             % output size (pixels)
%       meta.seg_flat_diam_px    % flat-to-flat hex size (pixels)
%       meta.focal_length_m      % focal length used for ideal focus (m)
%       meta.lambda_m            % wavelength (m)
%       meta.pixel_pitch_m       % pixel pitch in pupil sampling (m/pixel)
%
% Outputs
%   phi_full,  mask_full,  U_full      : assembled pupil (as before)
%   phi_stripped, mask_stripped, U_stripped :
%       same pupil after subtracting the exact global quadratic
%       matching meta.focal_length_m (best-focus pupil for FFT PSF).
%
% Notes
% - Requires hex_wavefront_random(..., seed, center_offset).

    % ---- pull geometry from meta ----
    img_res          = segments.meta.img_res;
    seg_flat_diam_px = segments.meta.seg_flat_diam_px;

    % ---- flat-top hex geometry / centers (pixels, origin at mosaic center) ----
    R  = seg_flat_diam_px / sqrt(3);    % center->vertex [px]
    aX = 1.5 * R;                        % axial->pixel x
    aY = sqrt(3) * R;                    % axial->pixel y

    axial = generate_axial_37();         % [q r] coordinates
    q = axial(:,1);  r = axial(:,2);

    Xc = aX * q;                         % mosaic-centered pixel coords
    Yc = aY * (r + q/2);

    % ---- place mosaic in an img_res x img_res frame, centered ----
    H = img_res;  W = img_res;
    cxG = (W+1)/2; cyG = (H+1)/2;
    XcG = cxG + Xc;                      % true (possibly fractional) centers
    YcG = cyG + Yc;

    % ---- outputs ----
    U_full    = zeros(H, W);
    phi_full  = zeros(H, W);
    mask_full = false(H, W);

    % ---- per-tile terms (fixed order) ----
    terms = {'piston','tilt_x','tilt_y','defocus'};

    % ---- loop over 37 segments ----
    for i = 1:37
        % per-segment coefficient means (radians); stds = 0 (deterministic)
        mu_sigma_i = [ segments.pistons(i), 0;
                       segments.tilt_x(i),  0;
                       segments.tilt_y(i),  0;
                       segments.defocus(i), 0 ];

        % true center & fractional offset (to remove seams)
        cx_true = XcG(i);  cy_true = YcG(i);
        cx_pix  = round(cx_true);  cy_pix  = round(cy_true);
        frac    = [cx_true - cx_pix,  cy_true - cy_pix];   % [dx, dy] pixels

        % generate cropped hex tile with sub-pixel center offset
        [U, phi, msk] = hex_wavefront_random(seg_flat_diam_px, terms, mu_sigma_i, [], frac);

        % paste by center alignment, FIRST-WINS to avoid border speckles
        h = size(msk,1);  w = size(msk,2);
        half_h = (h-1)/2;  half_w = (w-1)/2;

        xi0 = cx_pix - half_w;
        yi0 = cy_pix - half_h;
        rows = yi0 : yi0 + h - 1;
        cols = xi0 : xi0 + w - 1;

        % clip to the output frame
        r_keep = rows>=1 & rows<=H;
        c_keep = cols>=1 & cols<=W;
        if ~all(r_keep)
            rows = rows(r_keep);  msk = msk(r_keep,:);  U = U(r_keep,:);  phi = phi(r_keep,:);
        end
        if ~all(c_keep)
            cols = cols(c_keep);  msk = msk(:,c_keep);  U = U(:,c_keep);  phi = phi(:,c_keep);
        end
        if isempty(rows) || isempty(cols), continue; end

        already = mask_full(rows, cols);
        sub     = msk & ~already;

        Uf = U_full(rows, cols);   Uf(sub)  = U(sub);    U_full(rows, cols)  = Uf;
        Pf = phi_full(rows, cols); Pf(sub)  = phi(sub);  phi_full(rows, cols)= Pf;

        mask_full(rows, cols) = already | msk;
    end

    % ---- build the defocus STRIP term that matches the given focal length ----
    f0   = segments.meta.focal_length_m;
    lam  = segments.meta.lambda_m;
    pixm = segments.meta.pixel_pitch_m;

    R_m  = R * pixm;                                   % [m]
    if isinf(f0) || f0==0
        c0 = 0;
    else
        % Global quadratic coefficient for ideal focusing at f0
        c0 = -(pi/lam) * (R_m^2) * (1/f0);            % radians
    end

    % Global normalized coords (u = x/R, v = y/R) about the mosaic center
    [YY,XX] = ndgrid(1:H, 1:W);
    u = (XX - cxG) / R;
    v = (YY - cyG) / R;

    phi_lens = c0 * (u.^2 + v.^2);

    % Strip that exact quadratic (best-focus pupil for FFT)
    phi_stripped  = (phi_full - phi_lens) .* mask_full;
    mask_stripped = mask_full;
    U_stripped    = mask_full .* exp(1i*phi_stripped);
end
