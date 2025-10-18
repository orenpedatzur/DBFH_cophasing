% Define terms (include piston/tilts/defocus; add others if you want)
terms = {'piston','tilt_x','tilt_y','defocus'};

lambda   = 550e-9;     % m
pix      = 10e-6;      % m/pixel
seg_px   = 150;        % pixels
f0_m     = 120;        % mean focal length [m]
sigma_f  = 0.02;       % m (2 cm focal-length jitter around 120 m)
img_res = 2048 ; % square image edge size

% Tile physical radius (center->vertex)
R_m = (seg_px/sqrt(3)) * pix;

% --- DEFOCUS MEAN (radians): focusing curvature for f = f0_m ---
% Sign convention: converging focus => negative quadratic phase.
defocus_mean = -(pi/lambda) * (R_m^2) * (1/f0_m);

% --- DEFOCUS STD (radians): from sigma_f (meters) around f0_m ---
defocus_std  = defocus_std_from_sigma_f(seg_px, pix, lambda, sigma_f, f0_m);

% Base means/stds (radians) BEFORE per-tile tilt alignment.
mu_sigma = [ 0,            0.05;   % piston
             0,            0.02;   % tilt_x
             0,            0.02;   % tilt_y
             defocus_mean, defocus_std ];  % defocus: mean encodes f=120 m

seed = 1;

[U_full, phi_full, M_full, centers_uv, coeffs_all] = ...
    hex_aperture_wavefront(seg_px, terms, mu_sigma, f0_m, seed,img_res);

figure; imagesc(angle(U_full).*M_full);
axis image ij off; colormap(hsv); colorbar;
title('37-segment wavefront (phase, rad) — mean focus 120 m');


mask = segment_hex_mask_37([1 2 8 12 20 31 37], img_res, seg_px);
mask = double(mask);
mask(mask==0)=nan;

figure; imagesc(mask.*angle(U_full).*M_full);

% ---------- helper ----------
function c_std = defocus_std_from_sigma_f(seg_px, pixel_pitch_m, lambda_m, sigma_f_m, f0_m)
    if nargin < 5, f0_m = 120; end
    R_m  = (seg_px/sqrt(3)) * pixel_pitch_m;
    c_std = (pi/lambda_m) * (R_m^2 / (f0_m^2)) * sigma_f_m;
end
