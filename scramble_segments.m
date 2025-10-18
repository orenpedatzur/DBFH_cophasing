function segments_out = scramble_segments(segments, varargin)
%JITTER_SEGMENTS  Add random piston/tilt/defocus to each mirror segment.
%
% segments_out = jitter_segments(segments, Name, Value, ...)
%
% Inputs
%   segments : struct with fields (each length 37)
%       pistons, tilt_x, tilt_y, defocus
%     Optional: segments.meta with fields used for physical defocus jitter:
%       .seg_flat_diam_px, .pixel_pitch_m, .lambda_m, .focal_length_m
%
% Name-Value options (all radians unless noted)
%   'piston_std'   : default 0.05
%   'tilt_std_x'   : default 0.02
%   'tilt_std_y'   : default 0.02
%   'defocus_std'  : default 0.00   (ignored if 'sigma_f_m' provided)
%   'sigma_f_m'    : [] (meters). If provided, converts σ_f (m) -> defocus_std (rad)
%                    using segments.meta.{seg_flat_diam_px,pixel_pitch_m,lambda_m,focal_length_m}.
%   'seed'         : [] (set an integer for reproducibility)
%
% Output
%   segments_out : same as input but with random perturbations added.
%
% Notes
% - All jitters are added independently to each of the 37 segments.
% - If you pass 'sigma_f_m', we compute defocus_std = (π/λ) * (R_m^2 / f0^2) * σ_f,
%   where R_m = (seg_flat_diam_px/√3)*pixel_pitch_m and f0 = focal_length_m.

    % ---- parse options ----
    p = inputParser;
    p.addParameter('piston_std',  0.05, @(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
    p.addParameter('tilt_std_x',  0.02, @(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
    p.addParameter('tilt_std_y',  0.02, @(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
    p.addParameter('defocus_std', 0.00, @(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
    p.addParameter('sigma_f_m',   [],   @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isfinite(x)));
    p.addParameter('seed',        [],   @(x)isempty(x)||(isscalar(x)&&isfinite(x)));
    p.parse(varargin{:});
    opt = p.Results;

    if ~isempty(opt.seed)
        rng(opt.seed);
    end

    % ---- sizes & shape helpers ----
    N = numel(segments.pistons);
    assert(N==37, 'segments fields must have 37 elements.');
    col = @(v) reshape(v, [], 1);   % force column vectors

    % ---- resolve defocus_std (possibly from σ_f) ----
    defocus_std = opt.defocus_std;
    if ~isempty(opt.sigma_f_m)
        % need meta to convert meters -> radians
        need = {'seg_flat_diam_px','pixel_pitch_m','lambda_m','focal_length_m'};
        if ~isfield(segments,'meta') || any(~isfield(segments.meta, need))
            error(['sigma_f_m provided but segments.meta is missing required fields: ' ...
                   'seg_flat_diam_px, pixel_pitch_m, lambda_m, focal_length_m']);
        end
        seg_px   = segments.meta.seg_flat_diam_px;
        pix_m    = segments.meta.pixel_pitch_m;
        lambda_m = segments.meta.lambda_m;
        f0_m     = segments.meta.focal_length_m;

        R_m = (seg_px/sqrt(3)) * pix_m;
        defocus_std = (pi/lambda_m) * (R_m^2 / (f0_m^2)) * opt.sigma_f_m;
    end

    % ---- build randomized offsets ----
    dp  = opt.piston_std  * randn(N,1);
    dtx = opt.tilt_std_x  * randn(N,1);
    dty = opt.tilt_std_y  * randn(N,1);
    dd  = defocus_std     * randn(N,1);

    % ---- output struct ----
    segments_out = segments;  % copy meta and any other fields through
    segments_out.pistons = col(segments.pistons) + dp;
    segments_out.tilt_x  = col(segments.tilt_x)  + dtx;
    segments_out.tilt_y  = col(segments.tilt_y)  + dty;
    segments_out.defocus = col(segments.defocus) + dd;
end
