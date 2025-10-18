function out = dbh_prepare_from_four_centered(I1,I2,I12,I12p, dtilt1,dtilt2, Nfft, Rpx, pad_pixels)
%DBH_PREPARE_FROM_FOUR_CENTERED
% Build |A|, |B| and S = Â·B̂* from FOUR centered PSFs with known commanded tilts.
%
% Inputs (all HxW, CENTERED = fftshifted):
%   I1,I2      : intensities for tile 1 only, tile 2 only
%   I12        : intensity for tiles 1+2 together
%   I12p       : same as I12, but tile 2 has +pi/2 piston
%   dtilt1     : [tx ty] radians (linear-phase coeffs in u,v) applied to tile 1
%   dtilt2     : [tx ty] radians (linear-phase coeffs in u,v) applied to tile 2
%   Nfft       : FFT size used to compute these PSFs (defaults to size(I1,1))
%   Rpx        : hex radius (center->vertex) in pixels used in the pupil model
%   pad_pixels : scalar or [py px] zero-pad (default 32) to avoid circular wrap
%
% Output:
%   out.magA   = |Â|   (centered)
%   out.magB   = |B̂|  (centered)
%   out.S      = Â·B̂* (centered, complex)
%   out.C0, out.C90   (diagnostics)
%   out.I1_al, out.I2_al  (aligned component intensities, centered)
%   out.shifts_px = [dm1 dn1; dm2 dn2]  (applied image shifts in pixels)
%   out.Nfft, out.Rpx, out.pad

    if nargin < 7 || isempty(Nfft),       Nfft = size(I1,1); end
    if nargin < 8 || isempty(Rpx),        error('Provide Rpx (pixels).'); end
    if nargin < 9 || isempty(pad_pixels), pad_pixels = 0; end

    assert(isequal(size(I1),size(I2),size(I12),size(I12p)), 'Size mismatch among inputs.');

    % ---- tilt (rad in u,v) -> pixel shift in the CENTERED intensity image ----
    f = Nfft/(2*pi*Rpx);
    dm1 = dtilt1(1)*f;  dn1 = dtilt1(2)*f;
    dm2 = dtilt2(1)*f;  dn2 = dtilt2(2)*f;

    % ---- non-circular Fourier shift (centered spectra, with padding) ----
    I1_al = fourier_shift_centered_real_padded(I1,  dm1, dn1, pad_pixels);
    I2_al = fourier_shift_centered_real_padded(I2,  dm2, dn2, pad_pixels);

    % ---- DBH ingredients on a common, aligned grid (centered) ----
    C0  = I12  - I1_al - I2_al;   % ≈ 2*Re(Â·B̂*)
    C90 = I12p - I1_al - I2_al;   % ≈ 2*Im(Â·B̂*)
    S   = 0.5*(C0 + 1i*C90);      % Â·B̂*

    out = struct( ...
        'magA', sqrt(max(I1_al,0)), ...
        'magB', sqrt(max(I2_al,0)), ...
        'S',    S, ...
        'C0',   C0, ...
        'C90',  C90, ...
        'I1_al', I1_al, ...
        'I2_al', I2_al, ...
        'shifts_px', [dm1 dn1; dm2 dn2], ...
        'Nfft', Nfft, ...
        'Rpx',  Rpx, ...
        'pad',  pad_pixels );
end

function Yc = fourier_shift_centered_real_padded(Xc, dm, dn, pad_pixels)
% Shift a CENTERED (fftshifted) REAL image by (dm,dn) pixels via a frequency ramp,
% with zero padding to avoid circular wrap. Output remains CENTERED and real.

    if nargin < 4 || isempty(pad_pixels), pad_pixels = 32; end
    if isscalar(pad_pixels), pad_pixels = [pad_pixels pad_pixels]; end
    py = pad_pixels(1);  px = pad_pixels(2);

    % Quick exit for tiny shifts
    if max(abs([dm dn])) < 1e-12
        Yc = Xc;
        return;
    end

    [H,W] = size(Xc);

    % Centered pad
    Xc_big = padarray(Xc, [py px], 0, 'both');     % requires IPT; replace if needed
    [Hp,Wp] = size(Xc_big);

    % Centered spectrum of centered image
    Fc = fftshift( fft2( ifftshift(Xc_big) ) );

    % Centered frequency bins (−0.5..0.5 in steps of 1/Wp etc.)
    kx = -floor(Wp/2) :  ceil(Wp/2)-1;
    ky = -floor(Hp/2) :  ceil(Hp/2)-1;
    [KX,KY] = meshgrid(kx/Wp, ky/Hp);

    % Phase ramp: +dm → shift right, +dn → shift down
    ramp = exp(-1i*2*pi*( dm*KX + dn*KY ));

    % Apply ramp, inverse (centered) FFT
    Ybig = fftshift( ifft2( ifftshift(Fc .* ramp) ) );
    Ybig = real(Ybig);

    % Center crop back to HxW
    r0 = floor((Hp - H)/2) + 1;  c0 = floor((Wp - W)/2) + 1;
    Yc = Ybig(r0:r0+H-1, c0:c0+W-1);

    % Tiny numerical cleanup
    mx = max(1, max(abs(Yc(:))));
    tiny = 1e-12*mx;
    neg = (Yc < 0) & (Yc > -tiny);
    if any(neg,'all'), Yc(neg) = 0; end
end
