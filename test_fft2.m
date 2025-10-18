[E, I] = phase_fft2(phi, mask);
figure; imagesc(I); axis image ij off; colorbar
title('PSF (normalized intensity)');
