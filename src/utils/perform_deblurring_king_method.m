function [im_king, im_fs, tc, fc, fc_XYZ] = perform_deblurring_king_method(kspace, nufft_st, w, csm, gx, gy, gz, x, y, z, rotMatrixRCSToGCS, rotMatrixGCSToPCS, rotMatrixPCSToDCS, field_of_view_mm, DCS_offset, gamma, B0, dt)
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 12/27/2020, Last modified: 12/27/2020

%% Get data dimension
[Nk,Ni,Nc] = size(kspace);
[N1,N2,Nc] = size(csm);

%% Calculate the transformation matrix from RCS to DCS
%--------------------------------------------------------------------------
% From King 1999 MRM: xyz in the physical coordinates
% [x]   [a1 a2 a3][X]
% [y] = [a4 a5 a6][Y]
% [z]   [a7 a8 a9][Z]
%
% [gx]   [a1 a2 a3][GX]
% [gy] = [a4 a5 a6][GY]
% [gz]   [a7 a8 a9][GZ]
%
% For Siemens:
% [x]                                       [r]
% [y] = R_pcs2dcs * R_gcs2pcs * R_rcs2gcs * [c] + DCS_offset
% [z]                                       [s]
%
% xyz = A * rcs + A * A.' * DCS_offset
%     = A * (rcs + A.' * DCS_offset)
%     = A * XYZ
%
% [gx]                                       [gr]
% [gy] = R_pcs2dcs * R_gcs2pcs * R_rcs2gcs * [gc]
% [gz]                                       [gs]
%--------------------------------------------------------------------------
A = rotMatrixPCSToDCS * rotMatrixGCSToPCS * rotMatrixRCSToGCS;

a1 = A(1,1); a2 = A(1,2); a3 = A(1,3);
a4 = A(2,1); a5 = A(2,2); a6 = A(2,3);
a7 = A(3,1); a8 = A(3,2); a9 = A(3,3);

%% Calculate constants (F1 to F6)
F1 =  1 / 4 * (a1^2 + a4^2) * (a7^2 + a8^2) + a7^2 * (a2^2 + a5^2) - a7 * a8 * (a1 * a2 + a4 * a5);
F2 =  1 / 4 * (a2^2 + a5^2) * (a7^2 + a8^2) + a8^2 * (a1^2 + a4^2) - a7 * a8 * (a1 * a2 + a4 * a5);
F3 =  1 / 4 * (a3^2 + a6^2) * (a7^2 + a8^2) + a9^2 * (a1^2 + a2^2 + a4^2 + a5^2) - a7 * a9 * (a1 * a3 + a4 * a6) - a8 * a9 * (a2 * a3 + a5 * a6);
F4 =  1 / 2 * (a2 * a3 + a5 * a6) * (a7^2 - a8^2) + a8 * a9 * (2 * a1^2 + a2^2 + 2 * a4^2 + a5^2) - a7 * a8 * (a1 * a3 + a4 * a6) - a7 * a9 * (a1 * a2 + a4 * a5);
F5 =  1 / 2 * (a1 * a3 + a4 * a6) * (a8^2 - a7^2) + a7 * a9 * (a1^2 + 2 * a2^2 + a4^2 + 2 * a5^2) - a7 * a8 * (a2 * a3 + a5 * a6) - a8 * a9 * (a1 * a2 + a4 * a5);
F6 = -1 / 2 * (a1 * a2 + a4 * a5) * (a7^2 + a8^2) + a7 * a8 * (a1^2 + a2^2 + a4^2 + a5^2);

%% Calculate the spatial coordinates in RCS
xyz = cat(1, x(:).', y(:).', z(:).'); % 3 x N
XYZ = A.' * xyz;
X = reshape(XYZ(1,:), [N1 N2]);
Y = reshape(XYZ(2,:), [N1 N2]);
Z = reshape(XYZ(3,:), [N1 N2]);

%% Calculate the gradients in RCS
g_xyz = cat(1, gx(:).', gy(:).', gz(:).');
g_XYZ = A.' * g_xyz;
GX = reshape(g_XYZ(1,:), [Nk Ni]); % Nk x Ni [G/cm]
GY = reshape(g_XYZ(2,:), [Nk Ni]); % Nk x Ni [G/cm]
GZ = reshape(g_XYZ(3,:), [Nk Ni]); % Nk x Ni [G/cm]

%% Perform through-plane correction
start_time = tic; fprintf('Performing through-plane correction... ');
%--------------------------------------------------------------------------
% Calculate a scaled concomitant field time parameter tc(t) [sec]
%--------------------------------------------------------------------------
g0 = sqrt(GX.^2 + GY.^2); % Nk x Ni [G/cm]
gm = max(g0(:)); % [G/cm]
tc = 1 / gm^2 * cumsum(g0.^2 * dt); % 1 / [G/cm]^2 * [G/cm]^2 * [sec] => [sec]

%--------------------------------------------------------------------------
% Calculate the time-independent frequency offset fc
%--------------------------------------------------------------------------
XYZ_offset = A.' * DCS_offset;
% [rad/sec/T] / [2pi rad/cycle] * [m]^2 / [T] * ([G/cm] * [T/1e4G] * [1e2cm/m])^2 => [cycle/sec]
fc = gamma / (2 * pi) * (gm * 1e-2)^2 / (4 * B0) * F3 * XYZ_offset(3)^2; % [Hz]

%--------------------------------------------------------------------------
% Calculate the through-plane concomitant field phase
%--------------------------------------------------------------------------
phi_c = 2 * pi * fc * tc; % [rad/cycle] * [Hz] * [sec] => [rad]

%--------------------------------------------------------------------------
% Demodulation before gridding
%--------------------------------------------------------------------------
kspace_demod = bsxfun(@times, exp(-1j * phi_c), kspace);
fprintf('done! (%6.4f sec)\n', toc(start_time));

%% Perform NUFFT reconstruction
imc_nufft = complex(zeros(N1, N2, Nc, 'double'));
Nd = nufft_st.Nd;
for c = 1:Nc
    tic; fprintf('NUFFT reconstruction (c=%2d/%2d)... ', c, Nc);
    %----------------------------------------------------------------------
    % Apply the adjoint of 2D NUFFT
    %----------------------------------------------------------------------
    imc_nufft(:,:,c) = nufft_adj(reshape(double(kspace_demod(:,:,c)) .* w, [Nk*Ni 1]), nufft_st) / sqrt(prod(Nd));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end

%% IFFT to k-space (k-space <=> image-space)
kspace_demod_gridded = ifft2c(imc_nufft);

%% Perform in-plane correction (frequency-segmented deblurring)
L = 11; % number of frequency bins

%--------------------------------------------------------------------------
% Calculate the frequency offsets for concomitant fields
%--------------------------------------------------------------------------
fc_XYZ = gamma / (2 * pi) * (gm * 1e-2)^2 / (4 * B0) * (F1 * X.^2 + F2 * Y.^2 + F4 * Y .* Z + F5 * X .* Z + F6 * X .* Y); % [Hz]

%--------------------------------------------------------------------------
% Estimate the range of frequency offsets
% [rad/sec/T] / [2pi rad/cycle] * ([G/cm] * [T/1e4G] * [1e2cm/m])^2 / [T] * [m]^2 => [cycle/sec]
%--------------------------------------------------------------------------
D = max(field_of_view_mm) * 1e-3; % [mm] * [m/1e3mm] => [m]
index = (X.^2 + Y.^2 < (D / 2).^2);
df_max = max(fc_XYZ(index)); % [Hz]
interval = df_max / (L - 1); % [Hz]
df_range = (0:L-1).' * interval; % [Hz]
df_range = unique(df_range);
L = length(df_range);

%--------------------------------------------------------------------------
% Perform 2D interpolation on tc(t) [sec]
%--------------------------------------------------------------------------
om_scaled = nufft_st.om / (2 * pi);
KX_range = (-floor(N1/2):ceil(N1/2)-1).' * 1 / N1; % [-0.5,0.5]
KY_range = (-floor(N2/2):ceil(N2/2)-1).' * 1 / N2; % [-0.5,0.5]
[KX_grid,KY_grid,~] = ndgrid(KX_range, KY_range, 0);
F = scatteredInterpolant(om_scaled(:,1), om_scaled(:,2), double(tc(:)), 'linear', 'linear');
tc_grid = F(KX_grid, KY_grid);

%--------------------------------------------------------------------------
% Perform frequency-segmented deblurring
%--------------------------------------------------------------------------
im_fs = double(zeros(N1, N2, L, 'double'));
for idx1 = 1:L
    start_time = tic; fprintf('Performing frequency-segmented deblurring (%2d/%2d)... ', idx1, L);
    %----------------------------------------------------------------------
    % Demodulation after gridding
    %----------------------------------------------------------------------
    phi_c = 2 * pi * df_range(idx1) * tc_grid; % [rad/cycle] * [Hz] * [sec] => [rad]
    kspace_fs_gridded = bsxfun(@times, exp(-1j * phi_c), kspace_demod_gridded);

    %----------------------------------------------------------------------
    % FFT to image-space (k-space <=> image-space)
    %----------------------------------------------------------------------
    imc = crop(fft2c(kspace_fs_gridded), [N1 N2 Nc]);

    %----------------------------------------------------------------------
    % Perform optimal coil combination
    %----------------------------------------------------------------------
    im_fs(:,:,idx1) = sum(imc .* conj(csm), 3);
    fprintf('done! (%6.4f sec)\n', toc(start_time));
end

if L == 1
    im_king = im_fs;
    return;
end

%--------------------------------------------------------------------------
% Perform pixel-dependent interpolation (nearest neighbor)
%--------------------------------------------------------------------------
im_king = complex(zeros(N1, N2, 'double'));
for idx2 = 1:N2
    start_time = tic; fprintf('Performing pixel-dependent interpolation (:,%3d/%3d)... ', idx2, N2);
    for idx1 = 1:N1
        im_king(idx1,idx2) = interp1(df_range, reshape(im_fs(idx1,idx2,:), [L 1]), fc_XYZ(idx1,idx2), 'nearest', 'extrap');
    end
    fprintf('done! (%6.4f sec)\n', toc(start_time));
end

end