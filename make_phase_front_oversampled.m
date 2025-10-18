function [phi, mask, U, phi_os, mask_os, U_os] = make_phase_front_oversampled(segments, os)
%MAKE_PHASE_FRONT_OVERSAMPLED  Render pupil at os× resolution, anti-alias down.
%
% [phi, mask, U, phi_os, mask_os, U_os] = make_phase_front_oversampled(segments, os)
%
% Inputs
%   segments : your standard segments struct (base resolution)
%   os       : integer oversampling factor (e.g., 2 or 4)
%
% Outputs (base resolution)
%   phi, mask, U   : best-focus pupil at original img_res (downsampled)
% Also returns the oversampled best-focus pupil:
%   phi_os, mask_os, U_os : best-focus pupil at os*img_res
%
% Notes
% - Uses make_phase_front(segments_os) internally; expects it returns
%   [..., phi_stripped, mask_stripped, U_stripped] as the last three outputs.
% - Downsampling is anti-aliased via 'box' if available; otherwise uses
%   convolution with a box kernel and decimation.

    arguments
        segments struct
        os (1,1) {mustBeInteger, mustBePositive} = 2
    end

    % ---- Build an oversampled copy of the segment meta ----
    seg_os = segments;
    seg_os.meta.img_res            = os * segments.meta.img_res;
    seg_os.meta.seg_flat_diam_px   = os * segments.meta.seg_flat_diam_px;
    % Keep physical pupil size constant by shrinking pixel pitch:
    seg_os.meta.pixel_pitch_m      = segments.meta.pixel_pitch_m / os;

    % ---- Render the oversampled pupil (best-focus outputs) ----
    [~, ~, ~, phi_os, mask_os, U_os] = make_phase_front(seg_os);

    % ---- Anti-alias downsample back to base resolution ----
    H  = segments.meta.img_res;
    W  = segments.meta.img_res;

    % Complex-safe imresize
    if exist('imresize','file')
        try
            U = imresize(real(U_os), [H W], 'box') + 1i * imresize(imag(U_os), [H W], 'box');
            Maa = imresize(double(mask_os), [H W], 'box');       % anti-aliased mask in [0,1]
        catch
            % Older MATLAB: fall back if 'box' not supported
            [U, Maa] = box_downsample_fallback(U_os, mask_os, os, H, W);
        end
    else
        [U, Maa] = box_downsample_fallback(U_os, mask_os, os, H, W);
    end

    % Build hard mask from anti-aliased mask (threshold; tweak if desired)
    mask = Maa > 0.5;

    % Phase from complex field (masked)
    U  = mask .* U;              % zero out tiny edge leakage
    phi = angle(U);

    % ---- nested fallback (conv2+decimate) ----
    function [Uds, Mds] = box_downsample_fallback(Uhi, Mhi, s, Ht, Wt)
        k = ones(s,s) / s^2;
        % pad to avoid edge loss
        Uh = padarray(Uhi, [s s], 'replicate', 'both');
        Mh = padarray(double(Mhi), [s s], 'replicate', 'both');

        Uh = conv2(real(Uh), k, 'same') + 1i*conv2(imag(Uh), k, 'same');
        Mh = conv2(Mh, k, 'same');

        % center-preserving decimation indices
        % use imresize as reference mapping; here approximate:
        Uh = Uh( (1+s):s:end-s, (1+s):s:end-s );   % stride s
        Mh = Mh( (1+s):s:end-s, (1+s):s:end-s );

        % Crop or pad to exact [H W]
        Uds = imresize(real(Uh), [Ht Wt]) + 1i*imresize(imag(Uh), [Ht Wt]);
        Mds = imresize(Mh, [Ht Wt]);
    end
end
