function [piston_std, tilt_std, defocus_std] = stds_for_scramble(lambda, Rm, f, ...
                                           opd_nm, tilt_nrad, defocus_um)
k = 2*pi/lambda;
piston_std  = k * (opd_nm*1e-9);
tilt_std    = k * (2*tilt_nrad*1e-9) * Rm;       % same for x and y
defocus_std = (k*Rm^2/(2*f^2)) * (defocus_um*1e-6);
end
