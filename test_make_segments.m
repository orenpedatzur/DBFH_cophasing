close all;clear all;clc;

% Start from diffraction-limited segments
focal_distance = 120; % m
img_res = 2^10; 
seg_flat_diam_px = 140;
segments = make_segments(img_res, seg_flat_diam_px, focal_distance);

% Render to phase
[phi_0, mask,~,phi_0_stripped,~,~] = make_phase_front(segments);
imagesc(phi_0.*mask); axis image ij off;colorbar
title('diffraction limited phase front (rad)');

% Add random jitter (radians)
segments_noisy = scramble_segments(segments, ...
    'piston_std', 0.3, 'tilt_std_x', 0.3, 'tilt_std_y', 0.3,'defocus_std', 0.3,'seed', 1);

% Render to phase
[phi, ~,~,phi_stripped,~,~] = make_phase_front(segments_noisy);

% Fourier - Fraunhofer to focal plane
[E, I] = phase_fft2(phi_stripped, mask);
[E_0, I_0] = phase_fft2(phi_0_stripped, mask);

% display phases
subplot(3,3,3);
imagesc((phi-phi_0).*mask); axis image ij off;colorbar
title('Noisy exit-pupil phase - perfect phase (rad)');
subplot(3,3,1);
imagesc(phi_0.*mask); axis image ij off;colorbar
title('perfect phase (rad)');
subplot(3,3,2);
imagesc(phi.*mask); axis image ij off;colorbar
title('Noisy exit-pupil phase(rad)');

subplot(3,3,6);
imagesc((phi_stripped-phi_0_stripped).*mask); axis image ij off;colorbar
title('Noisy exit-pupil stripped phase - perfect phase (rad)');
subplot(3,3,5);
imagesc(phi_0_stripped.*mask); axis image ij off;colorbar
title('perfect stripped phase (rad)');
subplot(3,3,4);
imagesc(phi_stripped.*mask); axis image ij off;colorbar
title('Noisy exit-pupil stripped phase(rad)');

subplot(3,3,7); imagesc(log(I_0)); axis image ij off; colorbar
title('log PSF perfect wavefront(normalized intensity)');
subplot(3,3,8); imagesc(log(I)); axis image ij off; colorbar
title('log PSF noisy wavefront(normalized intensity)');

% check aliasing - maxStep needs to be < π
[maxStep, maxGrad] = report_phase_sampling(phi_stripped, mask);

% test unstacking
unstacked_segments = unstack_segment_tilts(segments_noisy);
[unstacked_phi, ~,~,unstacked_phi_stripped,~,~] = make_phase_front(unstacked_segments);

% unstacked @ Fourier - Fraunhofer
[unstacked_E, unstacked_I] = phase_fft2(unstacked_phi_stripped, mask);

% test pairing
pairs_1_segments = pairs_segments_tilt(unstacked_segments,1);
[pairs_1_phi, ~,~,pairs_1_phi_stripped,~,~] = make_phase_front(pairs_1_segments);

pairs_2_segments = pairs_segments_tilt(unstacked_segments,2);
[pairs_2_phi, ~,~,pairs_2_phi_stripped,~,~] = make_phase_front(pairs_2_segments);

% pairing @ Fourier - Fraunhofer
[pairs_1_E, pairs_1_I] = phase_fft2(pairs_1_phi_stripped, mask);
[pairs_2_E, pairs_2_I] = phase_fft2(pairs_2_phi_stripped, mask);



% display unstacking from a general position
figure;
subplot(3,4,1);
imagesc(phi.*mask); axis image ij off;colorbar
title('Noisy exit-pupil phase - perfect phase (rad)');
subplot(3,4,2);
imagesc(unstacked_phi.*mask); axis image ij off;colorbar
title('unstacked phase (rad)');
subplot(3,4,3);
imagesc(pairs_1_phi.*mask); axis image ij off;colorbar
title('pairs 1 phase (rad)');
subplot(3,4,4);
imagesc(pairs_2_phi.*mask); axis image ij off;colorbar
title('pairs 2 phase (rad)');

subplot(3,4,5);
imagesc(phi_stripped.*mask); axis image ij off;colorbar
title('Noisy exit-pupil stripped phase (rad)');
subplot(3,4,6);
imagesc(unstacked_phi_stripped.*mask); axis image ij off;colorbar
title('unstacked stripped phase (rad)');
subplot(3,4,7);
imagesc(pairs_1_phi_stripped.*mask); axis image ij off;colorbar
title('pairs 1 stripped phase (rad)');
subplot(3,4,8);
imagesc(pairs_2_phi_stripped.*mask); axis image ij off;colorbar
title('pairs 2 stripped phase (rad)');

subplot(3,4,9);
imagesc(log(I)); axis image ij off; colorbar
title('log PSF perfect wavefront(normalized intensity)');
subplot(3,4,10);
imagesc(log(unstacked_I)); axis image ij off; colorbar
title('log unstacked PSF wavefront(normalized intensity)');
subplot(3,4,11);
imagesc(log(pairs_1_I)); axis image ij off; colorbar
title('log pairs 1 PSF wavefront(normalized intensity)');
subplot(3,4,12);
imagesc(log(pairs_2_I)); axis image ij off; colorbar
title('log pairs 2 PSF wavefront(normalized intensity)');

sgtitle('test unstacking and pairing');



function [maxStep, maxGrad] = report_phase_sampling(phi, mask)
% maxStep: max absolute pixel-to-pixel difference (rad)
% maxGrad: max gradient magnitude via central differences (rad/pixel)
    phi(~mask) = NaN;
    dx = diff(phi,1,2); dy = diff(phi,1,1);
    maxStep = max([abs(dx(:)); abs(dy(:))], [], 'omitnan');

    % central-diff gradient
    Gx = (phi(:,[3:end end]) - phi(:,[1 1:end-2]))/2;
    Gy = (phi([3:end end],:) - phi([1 1:end-2],:))/2;
    maxGrad = max(sqrt(Gx(:).^2 + Gy(:).^2), [], 'omitnan');
end
