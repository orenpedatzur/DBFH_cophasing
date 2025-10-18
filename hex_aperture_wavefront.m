function [U_full, phi_full, mask_full, centers_uv, coeffs_all] = ...
    hex_aperture_wavefront(seg_flat_diam_px, terms, mu_sigma, focal_length, seed, img_res)
%HEX_APERTURE_WAVEFRONT
% Assemble a 37-segment wavefront by calling hex_wavefront_random().
% Per-tile mean piston  = value of all non-tilt MEAN terms at tile center.
% Per-tile mean tilts   = +gradient of those non-tilt MEAN terms (stitches smoothly).
% Center-aligned, first-wins pasting. Optional fixed output size via img_res.
% NEW: passes each tile's fractional center to hex_wavefront_random to avoid seams.

    if nargin < 5, seed = []; end
    if nargin < 6, img_res = []; end
    if ~isempty(seed), rng(seed); end %#ok<RNG>

    % ---------- mosaic geometry (flat-top) ----------
    R  = seg_flat_diam_px / sqrt(3);   % center->vertex [px]
    aX = 1.5 * R;                      % axial->pixel x
    aY = sqrt(3) * R;                  % axial->pixel y

    axial = generate_axial_37();       % [q r] (no drift / no y-skip)
    q = axial(:,1);  r = axial(:,2);

    % tile centers in pixels (origin at mosaic center)
    Xc = aX * q;
    Yc = aY * (r + q/2);

    % normalized centers used by hex_wavefront_random (u = x/R, v = y/R)
    Uc = Xc / R;  Vc = Yc / R;
    centers_uv = [Uc, Vc];

    % ---------- nominal per-tile size (matches flat-top hex) ----------
    % width = 2R (x), height = seg_flat_diam_px (y). Prefer odd sizes.
    nx0 = round(2*R);  ny0 = round(seg_flat_diam_px);
    if mod(nx0,2)==0, nx0 = nx0+1; end
    if mod(ny0,2)==0, ny0 = ny0+1; end
    half_w0 = (nx0-1)/2;  half_h0 = (ny0-1)/2;

    % ---------- choose global canvas ----------
    if ~isempty(img_res)
        % Fixed-size output, centered
        H = img_res;  W = img_res;
        cx_global = (W+1)/2;  cy_global = (H+1)/2;
        Xc_global = cx_global + Xc;   % true (possibly fractional) centers in frame
        Yc_global = cy_global + Yc;
    else
        % Auto-sized to fit all tiles (old behavior)
        pad = 1; % small safety margin
        xmin = floor(min(Xc) - half_w0) - pad;
        xmax = ceil( max(Xc) + half_w0) + pad;
        ymin = floor(min(Yc) - half_h0) - pad;
        ymax = ceil( max(Yc) + half_h0) + pad;
        W = xmax - xmin + 1;  H = ymax - ymin + 1;
        Xc_global = Xc - xmin + 1;
        Yc_global = Yc - ymin + 1;
    end

    % ---------- outputs ----------
    U_full     = zeros(H, W);
    phi_full   = zeros(H, W);
    mask_full  = false(H, W);
    coeffs_all = zeros(37, numel(terms));

    % term indices
    tnames  = string(lower(terms));
    id_pist = find(tnames=="piston",  1);
    id_tx   = find(tnames=="tilt_x",  1);
    id_ty   = find(tnames=="tilt_y",  1);

    base_mu = mu_sigma(:,1).';
    base_sg = mu_sigma(:,2).';

    % --------- per tile ---------
    for i = 1:37
        if ~isempty(seed), rng(seed + i); end %#ok<RNG>

        ui = Uc(i); vi = Vc(i);
        mu_i = base_mu;  sg_i = base_sg;

        % ---- mean phase & gradient at (ui,vi) from NON-tilt means ----
        phi_center = 0; gx = 0; gy = 0;
        for k = 1:numel(tnames)
            name = tnames(k);
            if name=="piston" || name=="tilt_x" || name=="tilt_y", continue; end
            c = mu_i(k);
            [dphi, dgx, dgy] = term_value_and_grad(name, c, ui, vi);
            phi_center = phi_center + dphi;
            gx = gx + dgx; gy = gy + dgy;
        end

        % ---- piston mean = value at center; tilt means = +gradient at center ----
        if ~isempty(id_pist), mu_i(id_pist) = mu_i(id_pist) + phi_center; end
        if ~isempty(id_tx),   mu_i(id_tx)   = mu_i(id_tx)   + gx;         end
        if ~isempty(id_ty),   mu_i(id_ty)   = mu_i(id_ty)   + gy;         end

        mu_sigma_i = [mu_i(:), sg_i(:)];

        % ---- true (possibly fractional) center in the output frame ----
        cx_true = Xc_global(i);
        cy_true = Yc_global(i);
        cx_pix  = round(cx_true);
        cy_pix  = round(cy_true);
        frac    = [cx_true - cx_pix,  cy_true - cy_pix];   % (dx,dy) in pixels

        % ---- generate tile with sub-pixel center offset ----
        % hex_wavefront_random signature: (..., seed, center_offset)
        [U, phi, msk, coeffs, ~] = hex_wavefront_random( ...
            seg_flat_diam_px, terms, mu_sigma_i, [], frac);
        coeffs_all(i,:) = coeffs;

        % ---- paste by CENTER alignment (first-wins) ----
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

        % first-wins pasting (prevents 1-px double-writes on borders)
        already = mask_full(rows, cols);
        sub     = msk & ~already;

        Uf = U_full(rows, cols);   Uf(sub) = U(sub);     U_full(rows, cols)  = Uf;
        Pf = phi_full(rows, cols); Pf(sub) = phi(sub);   phi_full(rows, cols)= Pf;

        mask_full(rows, cols) = already | msk;
    end
end

% ---------- helper: value & gradient of supported terms at (u,v) ----------
function [phi, gx, gy] = term_value_and_grad(name, c, u, v)
    switch name
        case "defocus"
            phi = c*(u*u + v*v);  gx = 2*c*u;  gy = 2*c*v;

        case "x2"
            phi = c*(u^2);        gx = 2*c*u;  gy = 0;
        case "xy"
            phi = c*(u*v);        gx = c*v;    gy = c*u;
        case "y2"
            phi = c*(v^2);        gx = 0;      gy = 2*c*v;

        % Cubics
        case "x3"
            phi = c*(u^3);        gx = 3*c*u^2;  gy = 0;
        case "x2y"
            phi = c*(u^2*v);      gx = 2*c*u*v;  gy = c*u^2;
        case "xy2"
            phi = c*(u*v^2);      gx = c*v^2;    gy = 2*c*u*v;
        case "y3"
            phi = c*(v^3);        gx = 0;        gy = 3*c*v^2;

        % Quartics
        case "x4"
            phi = c*(u^4);        gx = 4*c*u^3;  gy = 0;
        case "x3y"
            phi = c*(u^3*v);      gx = 3*c*u^2*v; gy = c*u^3;
        case "x2y2"
            phi = c*(u^2*v^2);    gx = 2*c*u*v^2; gy = 2*c*u^2*v;
        case "xy3"
            phi = c*(u*v^3);      gx = c*v^3;    gy = 3*c*u*v^2;
        case "y4"
            phi = c*(v^4);        gx = 0;        gy = 4*c*v^3;

        otherwise
            phi = 0; gx = 0; gy = 0;
    end
end
