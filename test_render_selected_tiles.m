% test_render_selected_tiles

close all;clear all;clc;

%% get phase in using 'make_phase_front'

% Start from diffraction-limited segments
focal_distance = 120; % m
img_res = 2^10; 
seg_flat_diam_px = 140;
segments = make_segments(img_res, seg_flat_diam_px, focal_distance);
unstacked_segments = unstack_segment_tilts(segments);
[unstacked_phi, mask,~,unstacked_phi_stripped,~,~] = make_phase_front(unstacked_segments);
[unstacked_E, unstacked_I] = phase_fft2(unstacked_phi_stripped, mask);

figure;
subplot(1,2,1);imagesc(unstacked_I)

%% get phase in using 'render_selected_tiles'

[phi37, M37, U37] = render_selected_tiles(unstacked_segments, [1:37]);
[~, I37]    = pupil_fft2(U37, img_res);

subplot(1,2,2);imagesc(I37)

%% both ways converge to the same result = good
