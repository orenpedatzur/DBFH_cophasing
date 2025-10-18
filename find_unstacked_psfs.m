function [x_psf, y_psf] = find_unstacked_psfs(unstacked_image)
%FIND_UNSTACKED_PSFS  Detect 37 spot peaks laid out on a hex-like grid.
%   [x_psf, y_psf] = FIND_UNSTACKED_PSFS(I) returns column vectors of
%   the peak coordinates (x = column index, y = row index) for image I.
%
%   The method:
%     1) light Gaussian smoothing to stabilize the peak,
%     2) regional-maximum detection with h-max suppression,
%     3) non-maximum suppression using a radius estimated from the data,
%     4) keep the top 37 peaks by intensity (enforcing one per spot).
%
%   Tips:
%     - Coordinates are in MATLAB pixel convention: (row=y, col=x).
%     - If your PSFs are blurrier/sharper, adjust 'sigma' and 'h_frac'.

% ---- parameters you can tweak if needed ----
sigma   = 5;     % Gaussian sigma (in pixels)
h_frac  = 0.10;    % h-max suppression as a fraction of dynamic range
max_keep = 37;     % expected number of spots
% --------------------------------------------

[n,~] = size(unstacked_image);

I = double(unstacked_image);
if ~ismatrix(I)
    error('Input must be a single-channel (grayscale) image.');
end

% Smooth & normalize
Is = imgaussfilt(I, sigma);
Is = Is - min(Is(:));
rngI = max(Is(:));
if rngI > 0, Is = Is ./ rngI; end

% find suppresion distance
f_Is = fft2(Is);
crop_ind = 3;
f_Is_max_search = f_Is(crop_ind:end/4,crop_ind:end/4);
[max_val,ind_max] = max(abs(f_Is_max_search),[],'all');
[max_ky,max_kx] = ind2sub(size(f_Is_max_search),ind_max);
max_kx = max_kx + crop_ind - 1;
max_ky = max_ky + crop_ind - 1;
max_k = sqrt(max_kx^2+max_ky^2);
diag_distance_px = n/cosd(30);
expected_pk2pk_dist = diag_distance_px/max_k;
supp_radius = expected_pk2pk_dist*0.75;

% Suppress shallow maxima, then find regional maxima
h = h_frac * (max(Is(:)) - min(Is(:)));
if h > 0
    Is2 = imhmax(Is, h);
else
    Is2 = Is;
end
BW = imregionalmax(Is2);             % candidate peak map

% Collect candidate coordinates & intensities
[r_all, c_all] = find(BW);
vals = Is2(sub2ind(size(Is2), r_all, c_all));

% If we somehow got too few candidates (over-suppression), relax h
if numel(vals) < max_keep
    BW = imregionalmax(Is);          % remove h-max, try again
    [r_all, c_all] = find(BW);
    vals = Is(sub2ind(size(Is), r_all, c_all));
end

% Sort candidates by brightness (desc)
[vals, idx] = sort(vals, 'descend');
r_all = r_all(idx);  c_all = c_all(idx);

% --- Estimate a reasonable suppression radius from data ---
% Use the K nearest neighbor distance among the brightest ~3*target
K = min(numel(vals), max_keep*3);
rK = r_all(1:K); cK = c_all(1:K);

% % For each of those, get nearest-neighbor distance (excluding itself)
% if K >= 2
%     D = pdist2([rK, cK], [rK, cK]);
%     D(1:K+1:end) = inf;                     % ignore self
%     nn = min(D, [], 2);                     % nearest neighbor dist
%     base_spacing = median(nn);              % robust grid spacing estimate
%     % suppression radius: smaller than spacing so we keep one per PSF
%     supp_radius = max(2, round(0.45 * base_spacing));
% else
%     % fallback if very few candidates (unlikely)
%     supp_radius = round( min(size(I))/20 );
% end

% --- Non-maximum suppression (greedy) ---
keep_r = zeros(max_keep, 1);
keep_c = zeros(max_keep, 1);
taken  = false(numel(vals), 1);
kept = 0;

for i = 1:numel(vals)
    if taken(i), continue; end
    kept = kept + 1;
    keep_r(kept) = r_all(i);
    keep_c(kept) = c_all(i);

    % suppress neighbors within radius
    dr = r_all - r_all(i);
    dc = c_all - c_all(i);
    taken = taken | (dr.^2 + dc.^2 <= supp_radius^2);

    if kept == max_keep
        break
    end
end

% If we still have fewer than 37, top up with the next best unsuppressed
if kept < max_keep
    remain = find(~taken);
    need = min(max_keep - kept, numel(remain));
    if need > 0
        extra = remain(1:need);
        keep_r(kept+1:kept+need) = r_all(extra);
        keep_c(kept+1:kept+need) = c_all(extra);
        kept = kept + need;
    end
end

% Return as (x,y) = (columns, rows), both column vectors
x_psf = keep_c(1:kept);
y_psf = keep_r(1:kept);

% Final safety: if more than 37 somehow, keep brightest 37 by sampling back
if numel(x_psf) > max_keep
    % compute intensities at kept points and trim
    v = Is2(sub2ind(size(Is2), y_psf, x_psf));
    [~, ii] = sort(v, 'descend');
    ii = ii(1:max_keep);
    x_psf = x_psf(ii);
    y_psf = y_psf(ii);
end

% ---------- ORDER PEAKS TO MATCH THE 1..37 HEX NUMBERING ----------
% Assumptions:
%   - index 1 is the center spot (closest to centroid)
%   - next 6 belong to the first ring, then 12, then 18
%   - within each ring, order starts at angle = 0 (pointing right/east)
%     and proceeds clockwise.

% Use the spots you've already selected
xk = nan([37,1]); yk = nan([37,1]);
xk(1:numel(x_psf)) = x_psf(:);  % columns (x)
yk(1:numel(x_psf)) = y_psf(:);  % rows (y)

% 1) robust center estimate (median helps if a point is a bit off)
xc = median(xk);
yc = median(yk);

% 2) radial distances and sort by radius
r = hypot(xk - xc, yk - yc);
[~, irad] = sort(r, 'ascend');
xk = xk(irad);  yk = yk(irad);  r = r(irad);

% 3) split into shells: 1 (center), 6, 12, 18
counts = [1, 6, 12, 18];
edges  = cumsum(counts);   % [1, 7, 19, 37]

order = zeros(37,1);   % indices (into current xk/yk) that realize 1..37
ptr = 1;

% --- center (1) ---
order(ptr) = 1;
ptr = ptr + 1;

% helper to sort a ring clockwise starting at +x
    function idx_ring_sorted = sort_ring(ix1, ix2)
        xr = xk(ix1:ix2);  yr = yk(ix1:ix2);
        % angle measured from +x, clockwise (note the minus sign on atan2)
        ang = atan2(yr - yc, xr - xc);
        % wrap to [0, 2pi)
        ang = mod(ang, 2*pi);

        % Start at angle closest to 0 (rightmost point), then clockwise
        [~, ii] = sort(ang, 'ascend');              % 0..2pi clockwise
        idx_ring_sorted = (ix1:ix2).';
        idx_ring_sorted = idx_ring_sorted(ii);
    end

% --- first ring (6) ---
order(ptr:ptr+counts(2)-1) = sort_ring(edges(1)+1, edges(2));
ptr = ptr + counts(2);

% --- second ring (12) ---
order(ptr:ptr+counts(3)-1) = sort_ring(edges(2)+1, edges(3));
ptr = ptr + counts(3);

% --- third ring (18) ---
order(ptr:ptr+counts(4)-1) = circshift(sort_ring(edges(3)+1, edges(4)),1);

% Apply the order
x_psf = xk(order);
y_psf = yk(order);
% ---------- END ORDERING ----------

end
