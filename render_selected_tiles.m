function [phi_best, mask_soft, U_best] = render_selected_tiles(segments, idx_on, edgeAA)
% Assemble only tiles in idx_on, create SOFT edge mask, strip the global quadratic,
% and return the complex pupil U_best (no re-binarization). edgeAA=1 px is good.

    if nargin < 3, edgeAA = 1; end

    img_res          = segments.meta.img_res;
    seg_flat_diam_px = segments.meta.seg_flat_diam_px;

    % mosaic geometry (keep your centers as-is; no integer forcing)
    R  = seg_flat_diam_px / sqrt(3);
    aX = 1.5*R;  aY = sqrt(3)*R;

    axial = generate_axial_37();  q=axial(:,1); r=axial(:,2);
    Xc = aX*q;  Yc = aY*(r + q/2);

    H = img_res; W = img_res;
    cxG = (W+1)/2;                       % keep fractional center if you want
    cyG = (H+1)/2;

    XcG = cxG + Xc;
    YcG = cyG + Yc;

    U_full = zeros(H,W);                 % complex pupil we'll FFT
    M_soft = zeros(H,W);                 % soft weights (for plotting/support)
    terms  = {'piston','tilt_x','tilt_y','defocus'};

    for ii = idx_on(:).'
        mu_sigma_i = [ segments.pistons(ii), 0;
                       segments.tilt_x(ii),  0;
                       segments.tilt_y(ii),  0;
                       segments.defocus(ii), 0 ];

        cx_true = XcG(ii);  cy_true = YcG(ii);
        cx_pix  = round(cx_true);  cy_pix = round(cy_true);
        frac    = [cx_true - cx_pix,  cy_true - cy_pix];

        [~, phi, msk] = hex_wavefront_random(seg_flat_diam_px, terms, mu_sigma_i, [], frac);

        % --- soft boundary: area-weighted pixel model (~1 px ramp) ---
        if edgeAA > 0
            Din  = bwdist(~msk);
            Dout = bwdist(msk);
            sd   = Din - Dout;                              % signed distance (px)
            msk_soft = max(0, min(1, 0.5 + sd/(2*edgeAA))); % linear ramp
        else
            msk_soft = double(msk);
        end

        Utile = msk_soft .* exp(1i*phi);    % complex tile with soft edge

        % paste (center alignment), first-wins to avoid seams
        [h, w] = size(msk_soft);
        half_h = (h-1)/2;  half_w = (w-1)/2;
        rows   = (cy_pix - half_h) : (cy_pix - half_h + h - 1);
        cols   = (cx_pix - half_w) : (cx_pix - half_w + w - 1);

        rk = rows>=1 & rows<=H;  ck = cols>=1 & cols<=W;
        rows = rows(rk); cols = cols(ck);
        if isempty(rows) || isempty(cols), continue; end
        Utile    = Utile(rk, ck);
        msk_soft = msk_soft(rk, ck);

        already = M_soft(rows, cols) > 0;
        take    = (~already) & (msk_soft > 0);

        Ublk = U_full(rows, cols);  Ublk(take) = Utile(take);
        U_full(rows, cols) = Ublk;

        Mblk = M_soft(rows, cols);  Mblk(take) = msk_soft(take);
        M_soft(rows, cols) = Mblk;
    end

    % strip the exact global quadratic at THE SAME center you used above
    f0   = segments.meta.focal_length_m;
    lam  = segments.meta.lambda_m;
    pixm = segments.meta.pixel_pitch_m;
    if isinf(f0) || f0==0
        c0 = 0;
    else
        c0 = -(pi/lam) * (R*pixm)^2 * (1/f0);   % phi_lens = c0*(u^2+v^2)
    end
    [YY,XX] = ndgrid(1:H, 1:W);
    u = (XX - cxG)/R;  v = (YY - cyG)/R;
    phi_lens = c0*(u.^2 + v.^2);

    U_best    = U_full .* exp(-1i*phi_lens);   % complex best-focus pupil
    mask_soft = M_soft;                        % keep soft weights
    phi_best  = angle(U_best);                 % phase if you want to visualize
    % figure;imagesc(real(U_full));axis equal ij
end
