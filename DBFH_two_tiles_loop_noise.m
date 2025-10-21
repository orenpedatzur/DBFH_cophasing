%% demo_dbfh_two_tiles.m
% Four captures for double-blind holography with tiles 1 & 2:
% 1) only tile 1,  2) only tile 2,
% 3) tiles 1+2 overlapped using NOMINAL tilts only,
% 4) tiles 1+2 overlapped + add piston pi/2 to tile 2.

clear all; clc;close all;
winsomnia(true);
parula_with_nan_white = [1 1 1; parula(256)];
font_size = 16;

%% Parameters
img_res = 2^11;
seg_px  = 140;     % flat-to-flat pixels per hex TODO: why only 140 works?
f0_m    = 120;     % focal length [m]
fft_res = img_res;    % PSF FFT size
seg_flat_diam_m = 6; % mirror size [m]
lambda = 550e-9;
k = 2*pi/lambda;
Rm = seg_flat_diam_m/sqrt(3);

%% Build diffraction-limited commands (uses your make_segments)
segments = make_segments(img_res, seg_px, f0_m);
% figure;scatter( segments.tilt_x, segments.tilt_y,'filled');axis ij equal
% [~, ~, U37] = render_selected_tiles(segments, [1:37]);
% [~, I37]    = pupil_fft2(U37, fft_res);
% 
% fig_unstacked = figure;imagesc((I37));  axis image  ij ; colorbar; colormap gray;
% title('scrambled and unstacked','Color','w');
%% add nontrivial phases to the segments
% pre-phasing tolerances:
sigma_opd_nm = 200; % nm
sigma_tilt_nrad = 50; % nrad
sigma_defocus_um = 5*50; % um

[piston_std, tilt_std, defocus_std] = stds_for_scramble(lambda, Rm, f0_m, ...
    sigma_opd_nm, sigma_tilt_nrad, sigma_defocus_um);

% scramble_segments expects tolerances in rads
segments = scramble_segments(segments,'piston_std', piston_std, 'tilt_std_x', tilt_std, 'tilt_std_y', tilt_std,'defocus_std', defocus_std,'seed', 4);

[phi37_scrambled, mask37_scrambled, ~] = render_selected_tiles(segments, [1:37]);
mask37_scrambled_nans = double(mask37_scrambled);
mask37_scrambled_nans(mask37_scrambled_nans==0) = nan;

figure; imagesc(phi37_scrambled.*mask37_scrambled_nans);  axis image ij off; colorbar; colormap(parula_with_nan_white);
set(gca,'FontSize',font_size);
exportgraphics(gcf,'figures\scrambled_phases.png');


%% unstack segments
unstack_ratio = 3.5e-5/16.4062 *f0_m*seg_px/img_res; % TODO: solve relationship between unstack_ratio and resolution & segment size
segments = unstack_segment_tilts(segments,unstack_ratio);
% figure;scatter( segments.tilt_x, segments.tilt_y,'filled');axis ij equal

[~, ~, U37] = render_selected_tiles(segments, [1:37]);
[~, I37]    = pupil_fft2(U37, fft_res);

% fig_unstacked = figure;imagesc((I37));  axis image  ij ; colorbar; colormap gray;
% title('scrambled and unstacked','Color','w');
% saveas(fig_unstacked,fullfile('figures','scrambled_and_unstacked'),'png');

test_numbering(fft_res,seg_px);
hold on;
imagesc(I37, 'AlphaData', 0.5);axis equal ij;colormap gray;
xlim([1 img_res]);
ylim([1 img_res]);

%% ---- Capture 1: only tile 1 ----
tile_1_ind = 3;
[phi1, M1, U1] = render_selected_tiles(segments, tile_1_ind);
[~, I1]    = pupil_fft2(U1, fft_res);

%% ---- Capture 2: only tile 2 ----
tile_2_ind = 8;
[phi2, M2, U2] = render_selected_tiles(segments, tile_2_ind);
[~, I2]    = pupil_fft2(U2, fft_res);

%% ---- Capture 3: tiles 1 & 2 overlapped using NOMINAL tilt only ----
% Align tile 2 spot to tile 1 spot in UNSTACKED tilt space, then re-stack.
[segments12,tilt_x,tilt_y] = overlap_by_nominal_tilt_unstacked(segments, tile_1_ind, tile_2_ind,unstack_ratio);
[phi12, M12, U12] = render_selected_tiles(segments12, [tile_1_ind tile_2_ind]);
[~, I12]     = pupil_fft2(U12, fft_res);

%% ---- Capture 4: tiles 1 & 2 overlapped + piston pi/2 on tile 2 ----
segments12p         = segments12;
segments12p.pistons(tile_2_ind) = segments12p.pistons(tile_2_ind) + pi/2;
[phi12p, M12p, U12p] = render_selected_tiles(segments12p, [tile_1_ind tile_2_ind]);
[~, I12p]      = pupil_fft2(U12p, fft_res);

%% Display (optional)
hex_grid = draw_hex_grid(segments);


figure;
tiledlayout(3,2,'Padding','tight','TileSpacing','compact');
nexttile; imagesc(log(I1));  axis image ij off; colormap parula; colorbar; title(['Tile ',num2str(tile_1_ind),' only']);set(gca,'FontSize',font_size);
nexttile; imagesc(log(I2));  axis image ij off; colormap parula; colorbar; title(['Tile ',num2str(tile_2_ind),' only']);set(gca,'FontSize',font_size);
nexttile; imagesc(log(I12)); axis image ij off; colormap parula; colorbar; title(['Tiles ',num2str(tile_1_ind),'+',num2str(tile_2_ind),' overlapped (nominal)']);set(gca,'FontSize',font_size);
nexttile; imagesc(log(I12p));axis image ij off; colormap parula; colorbar; title(['Tiles ',num2str(tile_1_ind),'+',num2str(tile_2_ind),' overlapped, piston \pi/2 on latter']);set(gca,'FontSize',font_size);
nexttile; imagesc(M1);hold on;imagesc(hex_grid,'AlphaData',0.3);axis image ij off; colormap parula; colorbar; title(['CS tile ',num2str(tile_1_ind)]);set(gca,'FontSize',font_size);
nexttile; imagesc(M2);hold on;imagesc(hex_grid,'AlphaData',0.3);axis image ij off; colormap parula; colorbar; title(['CS tile ',num2str(tile_2_ind)]);set(gca,'FontSize',font_size);
% exportgraphics(gcf,'figures\inputs_intensities_to_DBFH.png');


%% prepare data for DBH (two-tile case)
% After you generate the four captures (I1,I2,I12,I12p) and you KNOW the commanded tilts:
dtilt1 = [0 0];  % tile 1 unchanged
dtilt2 = [segments12.tilt_x(tile_2_ind) - segments.tilt_x(tile_2_ind), ...
    segments12.tilt_y(tile_2_ind) - segments.tilt_y(tile_2_ind)];
% dtilt2 = nominal_pair_dtilt(segments, tile_1_ind, tile_2_ind); % calc in the unstacked state

Rpx = segments.meta.seg_flat_diam_px / sqrt(3);


%% apply noise

peak_intensity_vec = 10.^[2:8];
sim_to_e = peak_intensity_vec ./ max(I12(:)); % scaling factor to move from simulation counts to realistic electron count-per-pixel


for ind_peak_instenisty = numel(peak_intensity_vec):-1:1

    this_sim_to_e = sim_to_e(ind_peak_instenisty);
    
    I1_noisy = add_jwst_noise(this_sim_to_e*I1, 1,'dark_rate_e',0,'read_noise_e', 0,'n_groups', 1,'bg_rate_e', 0);
    I2_noisy = add_jwst_noise(this_sim_to_e*I2, 1,'dark_rate_e',0,'read_noise_e', 0,'n_groups', 1,'bg_rate_e', 0);
    I12_noisy = add_jwst_noise(this_sim_to_e*I12, 1,'dark_rate_e',0,'read_noise_e', 0,'n_groups', 1,'bg_rate_e', 0);
    I12p_noisy = add_jwst_noise(this_sim_to_e*I12p, 1,'dark_rate_e',0,'read_noise_e', 0,'n_groups', 1,'bg_rate_e', 0);

    % Build |A|, |B|, S on the *centered* grid, then unshift to DFT layout
    prep = dbh_prepare_from_four_centered(I1_noisy, I2_noisy, I12_noisy, I12p_noisy, dtilt1, dtilt2, fft_res, Rpx);

    fprintf('max imag(I2_al) / max real(I2_al) = %.3g\n', ...
        max(abs(imag(prep.I2_al(:)))) / max(abs(real(prep.I2_al(:)))));
    fprintf('rel. RMS of |S|-|A||B|: %.3g\n', ...
        rms(abs(prep.S(:))-prep.magA(:).*prep.magB(:)) / mean(prep.magA(:).*prep.magB(:)));

    % --- DBH inputs in DFT layout (unshifted) ---
    magA = ifftshift(prep.magA);           % |Â|
    magB = ifftshift(prep.magB);           % |B̂|
    S    = ifftshift(prep.S);              % Â·B̂*  (complex)

    figure;
    subplot(1,3,1);imagesc(log(prep.magA));title('log(magA)');axis image ij off;
    subplot(1,3,2);imagesc(log(prep.magB));title('log(magB)');axis image ij off;
    subplot(1,3,3);imagesc(log(abs((prep.S))));title('log(abs((S)))');axis image ij off;

    % supports must be same size as magA/magB/S; center pad if needed
    suppA = center_padcrop(M1, size(magA), false);
    suppB = center_padcrop(M2, size(magA), false);

    % --- Build the linear system as operators (no huge explicit matrix) ---
    sys = dbh_prepare_system(magA, magB, S, suppA, suppB, 'lambda', 1e-3);

    % Real-stacked wrapper for lsqr (adapts complex A/A' to lsqr's API)
    Afun = @(x,tf) Afun_lsqr(x, tf, sys.Amul, sys.ATmul, sys.n, sys.m);

    % RHS and initial guess (real-stacked)
    bR = [real(sys.b); imag(sys.b)];
    x0 = zeros(2*sys.n,1);

    assert(numel(bR) == 2*sys.m, 'bR must be length 2*m');
    assert(numel(x0) == 2*sys.n, 'x0 must be length 2*n');

    % --- Solve ---
    tic;
    [xR,flag,relres,iter] = lsqr(Afun, bR, 1e-10, 1e5, [], [], x0);
    solve_time = toc;

    % Back to complex unknown on overlap K (phasor of B in frequency domain)
    zB = xR(1:sys.n) + 1i*xR(sys.n+1:end);
    zB = zB ./ max(abs(zB), 1e-12);   % (optional) unit-modulus projection

    % Recover zA on the overlap via S = A * conj(B): zA = (S / (|A||B|)) zB
    ZB = zeros(size(S));  ZB(sys.K) = zB;
    ZA = zeros(size(S));
    den = max(magA.*magB, 1e-15);
    ZA(sys.K) = (S(sys.K) ./ den(sys.K)) .* zB;

    %% --- Reconstruct pupil fields for the two tiles (DFT layout -> iFFT) ---
    Ahat_rec = magA .* ZA;      % Â_rec
    Bhat_rec = magB .* ZB;      % B̂_rec
    a_rec    = ifft2(Ahat_rec); % pupil field of tile 1 (complex)
    b_rec    = ifft2(Bhat_rec); % pupil field of tile 2 (complex)

    phi1_rec = angle(a_rec);    % recovered pupil phase (tile 1)
    phi2_rec = angle(b_rec);    % recovered pupil phase (tile 2)

    % display results (original vs. recovered, piston removed)
    M1 = logical(M1);
    M2 = logical(M2);

    % originals (your known pupil phases for tiles 1,2)
    phi1_0  = rm_piston_phi(phi1,  M1);
    % phi2_0  = rm_piston_phi(phi2,  M2);
    phi2_0  = rm_piston_phi(phi12.*M2,  M2);

    % recovered (from DBH)
    phi1r_0 = rm_piston_phi(phi1_rec, M1);
    phi2r_0 = rm_piston_phi(phi2_rec, M2);

    % differences & RMS (on their own masks)
    d1   = angle(exp(1i*(phi1r_0 - phi1_0))) .* M1;
    d2   = angle(exp(1i*(phi2r_0 - phi2_0))) .* M2;
    rms1(ind_peak_instenisty) = sqrt(mean(d1(M1).^2,'omitnan'));
    rms2(ind_peak_instenisty) = sqrt(mean(d2(M2).^2,'omitnan'));
    fprintf('RMS phase error (tile 1): %.3g rad\n', rms1(ind_peak_instenisty));
    fprintf('RMS phase error (tile 2): %.3g rad\n', rms2(ind_peak_instenisty));

    %% --- errors ---
    max_rms = max(rms1(ind_peak_instenisty),rms2(ind_peak_instenisty));
    piston_error_rms_m(ind_peak_instenisty) =  max_rms/sqrt(sum(M1(:))) / (2*pi) * (lambda)
    tilt_error_rms_rads(ind_peak_instenisty) =  2*max_rms*(lambda)/(2*pi)/sqrt(sum(M1(:)))/(seg_px/cosd(30)/2)
    defocus_error_rms_m(ind_peak_instenisty) = 4*sqrt(12)*(f0_m/(seg_flat_diam_m/2))^2*(max_rms*(lambda)/(2*pi))/sqrt(sum(M1(:)))

1
end

%% --- noise figure

gem = orderedcolors("gem");
max_rms = max([rms1;rms2],[],1);

figure;
subplot(4,1,1);loglog(SNR_vec,max_rms,'-o','LineWidth',2,'Color',gem(1,:));grid on;
xlabel('SNR (unitless)');ylabel('\phi_{px} RMSE [rad]');
subplot(4,1,2);loglog(SNR_vec,piston_error_rms_m,'-o','LineWidth',2,'Color',gem(2,:));grid on;
xlabel('SNR (unitless)');ylabel('piston RMSE [m]');
subplot(4,1,3);loglog(SNR_vec,tilt_error_rms_rads,'-o','LineWidth',2,'Color',gem(3,:));grid on;
xlabel('SNR (unitless)');ylabel('tilt x/y RMSE [rad]');
subplot(4,1,4);loglog(SNR_vec,defocus_error_rms_m,'-o','LineWidth',2,'Color',gem(4,:));grid on;
xlabel('SNR (unitless)');ylabel('defocus RMSE [m]');


%% --- Figure ---

[y, x] = find(M1);
xc_tile1 = mean(x);
yc_tile1 = mean(y);
[y, x] = find(M2);
xc_tile2 = mean(x);
yc_tile2 = mean(y);

nan_mask_tile1 = double(M1);
nan_mask_tile1(nan_mask_tile1==0) = nan;

nan_mask_tile2 = double(M2);
nan_mask_tile2(nan_mask_tile2==0) = nan;

cl1 = [-pi pi];
cl2 = 5e-2*[-pi pi];

font_size = 16;


figure;
tiledlayout(2,3, "TileSpacing", "compact", "Padding", "tight"); % Adjust TileSpacing and Padding

ax1 = nexttile; imagesc(phi1_0.*nan_mask_tile1);  axis image ij off; clim(cl1); colorbar; colormap(parula_with_nan_white);
title(['Tile  ',num2str(tile_1_ind),'  – original']);set(gca,'FontSize',font_size);


ax2 = nexttile; imagesc(phi1r_0.*nan_mask_tile1); axis image ij off; clim(cl1); colorbar;
title(['Tile ',num2str(tile_1_ind),' – recovered']);set(gca,'FontSize',font_size);

ax3 = nexttile; imagesc(d1.*nan_mask_tile1);      axis image ij off;clim(cl2);colorbar;
title(sprintf('Tile %d – diff (RMS=%.3g)',tile_1_ind,rms1));set(gca,'FontSize',font_size);

ax4 = nexttile; imagesc(phi2_0.*nan_mask_tile2);  axis image ij off; clim(cl1); colorbar;
title(['Tile  ',num2str(tile_2_ind),'  – original']);set(gca,'FontSize',font_size);

ax5 = nexttile; imagesc(phi2r_0.*nan_mask_tile2); axis image ij off; clim(cl1); colorbar;
title(['Tile  ',num2str(tile_2_ind),'  – recovered']);set(gca,'FontSize',font_size);

ax6 = nexttile; imagesc(d2.*nan_mask_tile2);      axis image ij off; clim(cl2);colorbar;
title(sprintf('Tile %d – diff (RMS=%.3g)',tile_2_ind,rms2));set(gca,'FontSize',font_size);


linkaxes([ax1,ax2,ax3]);
linkaxes([ax4,ax5,ax6]);

xlim(ax1,xc_tile1+0.8*seg_px*[-1,1]);
ylim(ax1,yc_tile1+0.8*seg_px*[-1,1]);

xlim(ax6,xc_tile2+0.8*seg_px*[-1,1]);
ylim(ax6,yc_tile2+0.8*seg_px*[-1,1]);

exportgraphics(gcf,'figures\DBFH_tiles_3_8.png');


%% --- end scripts ---
winsomnia(false);

%% ----------------- helpers (local functions) -----------------
function [segments_overlap,dtx,dty] = overlap_by_nominal_tilt_unstacked(unstacked_segments, ind1, ind2,unstack_ratio, fftN, signNom)
%OVERLAP_BY_NOMINAL_TILT_UNSTACKED  Use ONLY nominal tilts to overlap PSFs in UNSTACKED space.
%   segments_overlap = overlap_by_nominal_tilt_unstacked(unstacked_segments, ind1, ind2, fftN, signNom)
%
% Inputs
%   unstacked_segments : struct with residual (UNSTACKED) tilts
%   ind1, ind2         : segment indices (PSF of ind2 moves to ind1)
%   fftN               : optional FFT size for debug print (default = meta.img_res)
%   signNom            : optional sign for applying nominal delta; default = -1
%                        (i.e., tx2 <- tx2 - (tx_nom1 - tx_nom2))
%
% Output
%   segments_overlap   : same struct (still UNSTACKED); only ind2 tilts adjusted.

if nargin < 5 || isempty(fftN)
    fftN = unstacked_segments.meta.img_res;
end
if nargin < 6 || isempty(signNom)
    signNom = -1;   % your observed convention needs a minus
end

segments_overlap = unstacked_segments;
N = numel(unstacked_segments.tilt_x);
assert(ind1>=1 && ind1<=N && ind2>=1 && ind2<=N, 'ind1/ind2 out of range');

% --- geometry / nominal model ---
Rpx  = unstacked_segments.meta.seg_flat_diam_px / sqrt(3);
aX   = 1.5*Rpx;  aY = sqrt(3)*Rpx;
axial = generate_axial_37();
q = axial(:,1); r = axial(:,2);
Xc = aX*q; Yc = aY*(r + q/2);
u0 = Xc/Rpx; v0 = Yc/Rpx;

f0   = unstacked_segments.meta.focal_length_m;
lam  = unstacked_segments.meta.lambda_m;
pixm = unstacked_segments.meta.pixel_pitch_m;
if isinf(f0) || f0==0
    c0 = 0;
else
    Rm = Rpx * pixm;
    c0 = -(pi/lam) * (Rm^2) * (1/f0);   % phi = c0 (u^2+v^2)
end

tx_nom = unstack_ratio*2*c0*u0;   ty_nom = unstack_ratio*2*c0*v0;

% nominal delta: (1) - (2)
dtx = tx_nom(ind1) - tx_nom(ind2);
dty = ty_nom(ind1) - ty_nom(ind2);

% apply with chosen sign (default: minus)
tx = segments_overlap.tilt_x(:);
ty = segments_overlap.tilt_y(:);

% Predicted PSF shift (debug) before
kx_b = (tx/(2*pi*Rpx))*fftN;  ky_b = (ty/(2*pi*Rpx))*fftN;

tx(ind2) = tx(ind2) + signNom * dtx;
ty(ind2) = ty(ind2) + signNom * dty;

% Predicted PSF shift (debug) after
kx_a = (tx/(2*pi*Rpx))*fftN;  ky_a = (ty/(2*pi*Rpx))*fftN;

fprintf(['[overlap_nominal_unstacked] sign=%+d, Δt_nom=(%.4g,%.4g) rad  |  ' ...
    'PSF px BEFORE i1=(%.3f,%.3f) i2=(%.3f,%.3f)  ->  AFTER i1=(%.3f,%.3f) i2=(%.3f,%.3f)\n'], ...
    signNom, dtx, dty, ...
    kx_b(ind1), ky_b(ind1), kx_b(ind2), ky_b(ind2), ...
    kx_a(ind1),  ky_a(ind1),  kx_a(ind2),  ky_a(ind2));

segments_overlap.tilt_x = reshape(tx, size(unstacked_segments.tilt_x));
segments_overlap.tilt_y = reshape(ty, size(unstacked_segments.tilt_y));
end


function Y = center_padcrop(X, outSz, padval)
if nargin < 3, padval = 0; end
[h,w] = size(X);  H = outSz(1);  W = outSz(2);
Y = cast(padval*ones(H,W), 'like', X);
r0 = floor((H - h)/2);  c0 = floor((W - w)/2);
rs = max(1,1-r0) : min(h, H-r0);
cs = max(1,1-c0) : min(w, W-c0);
rd = (rs + r0);  cd = (cs + c0);
Y(rd, cd) = X(rs, cs);
end

function y = Afun_lsqr(x, transpFlag, Amul, ATmul, n, m)
% Real-stacked adapter for LSQR.
% - If transpFlag = 'notransp' (or 0): y = [Re;Im]( A * (Re+ i Im) )
% - If transpFlag = 'transp'   (or 1): y = [Re;Im]( A' * (Re+ i Im) )

if (ischar(transpFlag) && strcmpi(transpFlag,'transp')) || isequal(transpFlag,1)
    xr = x(1:m);  xi = x(m+1:end);
    yc = ATmul(xr + 1i*xi);       % complex m→n
    y  = [real(yc); imag(yc)];
else
    xr = x(1:n);  xi = x(n+1:end);
    yc = Amul(xr + 1i*xi);        % complex n→m
    y  = [real(yc); imag(yc)];
end
end

function phi0 = rm_piston_phi(phi, mask)
% Remove piston using circular mean over masked pixels.
% mask can be logical or numeric (0/1); outputs zero outside mask.
m = logical(mask);
if ~any(m,'all')
    phi0 = zeros(size(phi), 'like', phi);
    return
end
% circular mean phase on the mask
phibar = angle(mean(exp(1i*phi(m)), 'omitnan'));
% subtract and wrap to [-pi,pi), then zero outside mask
phi0 = angle(exp(1i*(phi - phibar))) .* m;
end

function dtilt2 = nominal_pair_dtilt(unstacked_segments, ind1, ind2)
% Return the NOMINAL tilt to move tile ind2's *unstacked* PSF onto tile ind1's,
% using only the nominal focusing model (no measured/residual terms).

% geometry
Rpx  = unstacked_segments.meta.seg_flat_diam_px / sqrt(3);
aX   = 1.5 * Rpx;
aY   = sqrt(3) * Rpx;
axial = generate_axial_37();
q = axial(:,1); r = axial(:,2);

Xc = aX * q;                   % pixels (mosaic centered)
Yc = aY * (r + q/2);
u  = Xc / Rpx;
v  = Yc / Rpx;

% nominal quadratic coefficient
f0   = unstacked_segments.meta.focal_length_m;
lam  = unstacked_segments.meta.lambda_m;
pixm = unstacked_segments.meta.pixel_pitch_m;
Rm   = Rpx * pixm;             % meters
if isinf(f0) || f0 == 0
    c0 = 0;
else
    c0 = -(pi/lam) * (Rm^2) * (1/f0);   % phi = c0*(u^2+v^2)
end

% delta tilt to move #2 onto #1 in the UNSTACKED plane (nominal only)
du = u(ind2) - u(ind1);
dv = v(ind2) - v(ind1);
dtilt2 = [ 2*c0*du,  2*c0*dv ];         % radians in the u,v basis
end

% % -test-
% function dtilt_nom = nominal_pair_dtilt(unstacked_segments, i1, i2)
%     Rpx  = unstacked_segments.meta.seg_flat_diam_px / sqrt(3);
%     aX   = 1.5*Rpx;  aY = sqrt(3)*Rpx;
%     axial = generate_axial_37();
%     q=axial(:,1); r=axial(:,2);
%     u = (aX*q)/Rpx;  v = (aY*(r+q/2))/Rpx;
%
%     f0   = unstacked_segments.meta.focal_length_m;
%     lam  = unstacked_segments.meta.lambda_m;
%     pixm = unstacked_segments.meta.pixel_pitch_m;
%     Rm   = Rpx * pixm;
%     c0   = (isinf(f0) || f0==0) ? 0 : -(pi/lam) * (Rm^2) * (1/f0);
%     % delta tilt of tile 2 to land on tile 1 in UNSTACKED space (nominal only)
%     dtilt_nom = [ 2*c0*(u(i1)-u(i2)), 2*c0*(v(i1)-v(i2)) ];
% end
