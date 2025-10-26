function noisy_e = add_jwst_noise(img_e, t_exp, opts)
% img_e : expected signal electrons (per pixel) for the full exposure
% t_exp : exposure time [s]
% opts  : struct with fields (set per instrument):
%   .dark_rate_e  (e-/s/pix)   e.g., 0.01 for NIRCam, 0.2 for MIRI
%   .read_noise_e (e- RMS/read) e.g., 12 for NIRCam, 14 for MIRI
%   .n_groups     (# UTR groups; default 1 = CDS)
%   .bg_rate_e    (e-/s/pix) optional astrophysical background
%
% Returns noisy electrons per pixel (same size as img_e).

    arguments
        img_e {mustBeNonnegative}
        t_exp (1,1) {mustBePositive}
        opts.dark_rate_e (1,1) {mustBeNonnegative} = 0.01      % NIRCam default
        opts.read_noise_e (1,1) {mustBeNonnegative} = 12        % NIRCam default
        opts.n_groups (1,1) {mustBeInteger, mustBePositive} = 1
        opts.bg_rate_e (1,1) {mustBeNonnegative} = 0
        opts.seed (1,1) = 1;
    end

    % expectation per pixel (electrons)
    mu = img_e + (opts.bg_rate_e + opts.dark_rate_e) * t_exp;

    % Poisson shot noise
    % (Protect against extremely large values; cast to double afterwards)
    rng(opts.seed);gpurng(opts.seed)
    noisy_poiss = poissrnd(mu);

    % Effective read noise for UTR (simple 1/sqrt(N) rule of thumb)
    sigma_r_eff = opts.read_noise_e / sqrt(opts.n_groups);

    % Add Gaussian read noise
    noisy_e = round(noisy_poiss + sigma_r_eff .* randn(size(img_e)));
    noisy_e(noisy_e<0) = 0;
end
