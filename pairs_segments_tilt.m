function tilted_segments = pairs_segments_tilt(segments, pairing_scheme)
%PAIRS_SEGMENTS_TILT  Pair segments in UNSTACKED tilt space and move each pair
% to their midpoint using only tilts, then re-stack.
%
% tilted_segments = pairs_segments_tilt(segments, pairing_scheme)
%
% pairing_scheme:
%   1 -> (2,3), (4,5), ..., (36,37)            % even with next odd
%   2 -> (3,4), (5,6), ..., (35,36)            % odd with next even; don't tilt 37
% Segment 1 is excluded in both schemes. Pistons/defocus/meta unchanged.
%
% Requires shared generate_axial_37() on path.

    if nargin < 2 || isempty(pairing_scheme), pairing_scheme = 1; end
    assert(ismember(pairing_scheme,[1 2]), 'pairing_scheme must be 1 or 2.');

    % ---- checks ----
    need = {'img_res','seg_flat_diam_px','focal_length_m','lambda_m','pixel_pitch_m'};
    assert(isfield(segments,'meta') && all(isfield(segments.meta, need)), ...
        'segments.meta must contain: %s', strjoin(need, ', '));

    N = numel(segments.tilt_x);
    assert(N == 37, 'segments.* fields must have 37 elements.');

    % preserve original shape (row/col)
    shp  = size(segments.tilt_x);
    tx0  = segments.tilt_x(:);
    ty0  = segments.tilt_y(:);

    % ---- geometry (normalized centers u0,v0) ----
    Rpx  = segments.meta.seg_flat_diam_px / sqrt(3);
    aX   = 1.5 * Rpx;  aY = sqrt(3) * Rpx;
    axial = generate_axial_37();
    q = axial(:,1);  r = axial(:,2);
    Xc = aX * q;  Yc = aY * (r + q/2);
    u0 = Xc / Rpx;  v0 = Yc / Rpx;

    % ---- nominal focusing tilts from focal length ----
    f0   = segments.meta.focal_length_m;
    lam  = segments.meta.lambda_m;
    pixm = segments.meta.pixel_pitch_m;

    if isinf(f0) || f0 == 0
        c0 = 0;
    else
        Rm = Rpx * pixm;                         % meters
        c0 = -(pi/lam) * (Rm^2) * (1/f0);        % radians (phi = c0*(u^2+v^2))
    end
    tx_nom = 2*c0*u0;                            % ∂phi/∂u at center
    ty_nom = 2*c0*v0;                            % ∂phi/∂v at center

    % ---- unstack to residual tilt space ----
    tx_un = tx0 - tx_nom;
    ty_un = ty0 - ty_nom;

    % ---- choose pairing indices ----
    switch pairing_scheme
        case 1  % even -> next odd : (2,3),(4,5),...,(36,37)
            idxA = (2:2:36).';          % even
            idxB = idxA + 1;            % next odd
        case 2  % odd  -> next even: (3,4),(5,6),...,(35,36); don't tilt 37
            idxA = (3:2:35).';          % odd (exclude 1 and 37)
            idxB = idxA + 1;            % next even
            % 2 and 37 remain unpaired/unchanged
    end

    % ---- pairwise midpoint in unstacked space ----
    tx_mid = (tx_un(idxA) + tx_un(idxB)) / 2;
    ty_mid = (ty_un(idxA) + ty_un(idxB)) / 2;

    tx_un(idxA) = tx_mid;  tx_un(idxB) = tx_mid;
    ty_un(idxA) = ty_mid;  ty_un(idxB) = ty_mid;

    % ---- re-stack (back to commanded tilt space) ----
    tx_new = tx_un + tx_nom;
    ty_new = ty_un + ty_nom;

    % ---- build output ----
    tilted_segments        = segments;
    tilted_segments.tilt_x = reshape(tx_new, shp);
    tilted_segments.tilt_y = reshape(ty_new, shp);

    % optional: warn if nothing changed
    tol = 1e-12;
    if all(abs(tilted_segments.tilt_x(:) - tx0) < tol) && ...
       all(abs(tilted_segments.tilt_y(:) - ty0) < tol)
        warning(['pairs_segments_tilt: no change to tilts. Likely all paired ' ...
                 'segments already had identical UNSTACKED tilts (e.g., stacked state).']);
    end
end
