function [E, I] = pupil_fft2(U, fft_res)
    if nargin < 2 || isempty(fft_res), fft_res = size(U,1); end
    % center-pad to fft_res
    padH = max(0, fft_res - size(U,1));
    padW = max(0, fft_res - size(U,2));
    Upad = padarray(U, [floor(padH/2), floor(padW/2)], 0, 'pre');
    Upad = padarray(Upad, [ceil(padH/2),  ceil(padW/2)],  0, 'post');
    E = fftshift(fft2(ifftshift(Upad),fft_res,fft_res));
    I = abs(E).^2;
end
