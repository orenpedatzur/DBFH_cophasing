function [E, I, U] = phase_fft2(phi, mask, fft_res)
%PHASE_FFT2  2-D FFT of a pupil phase map.
%   [E, I, U] = phase_fft2(phi)
%   [E, I, U] = phase_fft2(phi, mask)
%   [E, I, U] = phase_fft2(phi, mask, fft_res)
%
% Inputs
%   phi      : phase (radians), size HxW. Center assumed at (H+1)/2,(W+1)/2.
%   mask     : optional logical mask HxW (default = true(size(phi))).
%   fft_res  : optional scalar or [Hout Wout] for FFT size (zero-padding).
%              Default = size(phi).
%
% Outputs
%   E : complex field in focal plane, centered (fftshifted)
%   I : intensity |E|^2 normalized to peak = 1
%   U : complex pupil field used: U = mask .* exp(1i*phi)
%
% Notes
% - Uses ifftshift/fftshift so the pupil center in phi/mask maps to
%   the array center in E/I.
% - Normalizes E by sqrt(sum(mask(:))) so total energy is consistent,
%   then normalizes I to max(I)=1 for easy visualization.

    if nargin < 2 || isempty(mask)
        mask = true(size(phi));
    end
    if nargin < 3
        fft_res = [];
    end

    [H, W] = size(phi);

    % Output size
    if isempty(fft_res)
        Hout = H; Wout = W;
    elseif isscalar(fft_res)
        Hout = fft_res; Wout = fft_res;
    else
        Hout = fft_res(1);
        Wout = fft_res(2);
    end

    % Complex pupil field
    U = double(mask) .* exp(1i*phi);

    % Fraunhofer diffraction (frequency-domain field), centered
    E = fftshift(fft2(ifftshift(U), Hout, Wout));

    % Energy normalization (optional but standard)
    denom = sqrt(max(1, sum(mask(:))));
    E = E / denom;

    % Intensity, normalized to peak = 1
    I = abs(E).^2;
    I = I / max(I(:) + eps);
end
