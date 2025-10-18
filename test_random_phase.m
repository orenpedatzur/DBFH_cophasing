%% Parameters
lambda = 550e-9;          % wavelength [m]
f0 = 120;                 % nominal focus [m]
pixel_pitch = 10e-6;      % [m/pixel]
seg_px = 150;             % flat-to-flat size [pixels]

%% Random physical parameters
tilt_mean_deg = 0; tilt_std_deg = 0.1;       % random tilt ±0.01°
defocus_mean_m = 0; defocus_std_m = 0.2;     % random Δf ±2 cm

% Draw random angles/defocus
theta_x = deg2rad(tilt_mean_deg + tilt_std_deg*randn);
theta_y = deg2rad(tilt_mean_deg + tilt_std_deg*randn);
delta_f = defocus_mean_m + defocus_std_m*randn;

%% Convert to equivalent phase coefficients for hex_wavefront_random
% Define coordinate scale (center-to-vertex radius)
R = seg_px / sqrt(3);
half_size_m = R * pixel_pitch;    % physical half-width of segment [m]

% Estimate phase slopes (radians of phase per normalized coordinate unit)
% Normalized u = x / R, so full physical slope multiplied by R to normalize
tilt_coeff_x = (2*pi/lambda) * half_size_m * sin(theta_x);
tilt_coeff_y = (2*pi/lambda) * half_size_m * sin(theta_y);

% Defocus coefficient: curvature (π/(λ Δf)) * (half_size_m)^2
defocus_coeff = (pi/(lambda*delta_f)) * (half_size_m)^2;

% Random piston term (arbitrary)
piston_std = 0.2;  % ±0.1 rad typical
piston_coeff = 2*pi*randn*piston_std; 

%% Prepare term names and coefficient stats
terms   = {'piston','tilt_x','tilt_y','defocus'};
mu_sigma = [piston_coeff, 0;    % no randomness in this run
             tilt_coeff_x, 0;
             tilt_coeff_y, 0;
             defocus_coeff, 0];

%% Call your existing function
[U, phi, mask, coeffs, terms_out] = hex_wavefront_random(seg_px, terms, mu_sigma);

%% Display
figure;
imagesc(phi.*mask);
axis image ij off;
colormap(hsv); colorbar;
title(sprintf('Tilt=(%.4f°, %.4f°), Δf = %.3f m', ...
    rad2deg(theta_x), rad2deg(theta_y), delta_f));
