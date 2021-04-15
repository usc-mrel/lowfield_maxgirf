function [im_echo,B0map_fit_mask,x_echo,y_echo,z_echo,TE,csm] = cartesian_B0map_recon(gre_noise_fullpaths, gre_data_fullpaths, gre_dat_fullpaths, user_opts_cartesian)
% Inputs
%   gre_noise_full_paths    cell structure containing full paths to noise only h5 files
%   gre_data_full_paths     cell structure containing full paths to mutlti-echo h5 files
%   gre_dat_full_paths      cell structure containing full paths to Siemens dat files


%% Parse the user options
output_directory = user_opts_cartesian.output_directory;
zpad_factor      = user_opts_cartesian.zpad_factor;
main_orientation = user_opts_cartesian.main_orientation;
lambda           = user_opts_cartesian.lambda;
maxiter_tgv      = user_opts_cartesian.maxiter_tgv;

%% Read the first ismrmrd file to get header information
start_time = tic;
tic; fprintf('Reading an ismrmrd file: %s... ', gre_data_fullpaths{1});
if exist(gre_data_fullpaths{1}, 'file')
    dset = ismrmrd.Dataset(gre_data_fullpaths{1}, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , gre_data_fullpaths{1});
end

%% Parse the ISMRMRD header
raw_data = dset.readAcquisition(); % read all the acquisitions

%% Get imaging parameters from the XML header
header = ismrmrd.xml.deserialize(dset.readxml);

%% Set the dimension of the data
N1 = zpad_factor * header.encoding.reconSpace.matrixSize.x; % number of voxels in image-space (RO,row)
N2 = zpad_factor * header.encoding.reconSpace.matrixSize.y; % number of voxels in image-space (PE,col)
N3 = double(max(raw_data.head.idx.slice)) + 1; % number of slices
Ne = length(gre_data_fullpaths); % number of echoes
N  = N1 * N2; % total number of voxels in image-space

%% Reconstruct multi-echo GRE data
im_echo = complex(zeros(N1, N2, N3, Ne, 'double'));
TE = zeros(Ne, 1, 'double'); % [sec]

for idx = 1:Ne
    %----------------------------------------------------------------------
    % Perform Cartesian reconstruction
    %----------------------------------------------------------------------
    if idx == 1
        [im,header,x,y,z,imc,csm] = cartesian_recon(gre_noise_fullpaths{idx}, gre_data_fullpaths{idx}, gre_dat_fullpaths{idx}, user_opts_cartesian);
        csm = squeeze(csm); % N1 x N2 x N3 x Nc
    else
        [~,header,x,y,z,imc] = cartesian_recon(gre_noise_fullpaths{idx}, gre_data_fullpaths{idx}, gre_dat_fullpaths{idx}, user_opts_cartesian);
        imc = squeeze(imc); % N1 x N2 x N3 x Nc
        im = sum(conj(csm) .* imc, ndims(imc)) ./ sum(abs(csm).^2, ndims(imc));
    end

    %----------------------------------------------------------------------
    % Save the reconstructed image
    %----------------------------------------------------------------------
    im_echo(:,:,:,idx) = reshape(im, [N1 N2 N3]);
    TE(idx) = header.sequenceParameters.TE * 1e-3; % [msec] * [sec/1e3msec] => [sec]
end

%% Get spatial coordinates in DCS [m]
x_echo = reshape(x, [N1 N2 N3]); % N1 x N2 x N3
y_echo = reshape(y, [N1 N2 N3]); % N1 x N2 x N3
z_echo = reshape(z, [N1 N2 N3]); % N1 x N2 x N3

%% Calculate a raw static off-resonance map [Hz]
B0map_raw = zeros(N1, N2, N3, 'double'); % [Hz]

for idx3 = 1:N3
    tic; fprintf('(%2d/%2d): Calculating a static off-resonance map... ', idx3, N3);
    nr_voxels = N1 * N2;
    for idx = 1:nr_voxels
        %------------------------------------------------------------------
        % Perform least-squares fitting per voxel
        % m = a * TE + b, m in [rad], a in [rad/sec], TE in [sec]
        % m(1) = a * TE(1) + b    [m(1)] = [TE(1) 1][a]
        % m(2) = a * TE(2) + b => [m(2)] = [TE(2) 1][b] => m = Av
        % m(N) = a * TE(N) + b    [m(N)] = [TE(N) 1]
        %------------------------------------------------------------------
        [idx1,idx2] = ind2sub([N1 N2], idx);
        m = unwrap(reshape(angle(im_echo(idx1,idx2,idx3,:)), [Ne 1])); % [rad]
        A = [TE ones(Ne,1)]; % [sec]
        v = A \ m; % [rad/sec]
        B0map_raw(idx1,idx2,idx3) = v(1) / (2 * pi); % [rad/sec] * [cycle/2pi rad] => [Hz]
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end
    
%% Calculate a mask defining the image and noise regions using the iterative intermean algorithm
%--------------------------------------------------------------------------
% Calculate the maximum magnitude of images from all echo points
%--------------------------------------------------------------------------
im_max = max(abs(im_echo), [], 4);

%--------------------------------------------------------------------------
% Calculate a mask
%--------------------------------------------------------------------------
mask = false(N1, N2, N3);
for idx = 1:N3
    level = isodata(im_max(:,:,idx), 'log');
    mask(:,:,idx) = (im_max(:,:,idx) > level);
end

%% Display results
if main_orientation == 0 % sagittal plane
    reorient = @(x) x;
elseif main_orientation == 1 % coronal plane
    reorient = @(x) x;
elseif main_orientation == 2 % transverse plane
    reorient = @(x) rot90(x, -1);
end

figure('Color', 'w');
montage(reorient(im_max), 'DisplayRange', []);
export_fig(fullfile(output_directory, 'im_max'), '-r400', '-tif');

figure('Color', 'w');
montage(reorient(mask));
export_fig(fullfile(output_directory, 'mask'), '-r400', '-tif');

figure('Color', 'w');
montage(reorient(B0map_raw), 'DisplayRange', []);
colormap(hot(256)); colorbar;
caxis([-60 60]);
title('Raw static off-resonance map [Hz]')
export_fig(fullfile(output_directory, 'B0map'), '-r400', '-tif');

figure('Color', 'w');
montage(reorient(B0map_raw .* mask), 'DisplayRange', []);
colormap(hot(256)); colorbar;
caxis([-60 60]);
title('Masked raw static off-resonance map [Hz]')
export_fig(fullfile(output_directory, 'B0map_mask'), '-r400', '-tif');

%% Calculate a smooth spherical harmonics approximation to the off-resonance map
B0map_fit = zeros(N1, N2, N3, 'double');

for idx3 = 1:N3
    %% Set spatial coordinates of the current slice
    x = reshape(x_echo(:,:,idx3), [N 1]);
    y = reshape(y_echo(:,:,idx3), [N 1]);
    z = reshape(z_echo(:,:,idx3), [N 1]);

    %% Calculate real-valued spherical harmonic basis functions
    tic; fprintf('(%2d/%2d): Calculating real-valued spherical harmonic basis functions... ', idx3, N3);
    %----------------------------------------------------------------------
    % Calculate the number of basis functions
    %----------------------------------------------------------------------
    max_order = 5;
    Ns = 0; % number of basis functions
    for order = 0:max_order
        Ns = Ns + 2 * order + 1;
    end
    A = zeros(N, Ns, 'double');

    %----------------------------------------------------------------------
    % Order = 0 (1)
    %----------------------------------------------------------------------
    if (max_order >= 0)
        A(:,1) = 1;
    end

    %----------------------------------------------------------------------
    % Order = 1 (3)
    %----------------------------------------------------------------------
    if (max_order >= 1)
        A(:,2) = x;
        A(:,3) = y;
        A(:,4) = z;
    end

    %----------------------------------------------------------------------
    % Order = 2 (5)
    %----------------------------------------------------------------------
    if (max_order >= 2)
        A(:,5) = x .* y;
        A(:,6) = z .* y;
        A(:,7) = 2 * z.^2 - (x.^2 + y.^2);
        A(:,8) = x .* z;
        A(:,9) = x.^2 - y.^2;
    end

    %----------------------------------------------------------------------
    % Order = 3 (7)
    %----------------------------------------------------------------------
    if (max_order >= 3)
        A(:,10) = 3 * y .* x.^2 - y.^3;
        A(:,11) = x .* y .* z;
        A(:,12) = 5 * y .* z.^2 - y .* (x.^2 + y.^2 + z.^2);
        A(:,13) = 2 * z.^3 - 3 * z .* (x.^2 + y.^2);
        A(:,14) = 5 * x .* z.^2 - x .* (x.^2 + y.^2 + z.^2);
        A(:,15) = z .* (x.^2 - y.^2);
        A(:,16) = x.^3 - 3 * x .* y.^2;
    end

    %----------------------------------------------------------------------
    % Order = 4 (9)
    %----------------------------------------------------------------------
    if (max_order >= 4)
        r2 = x.^2 + y.^2 + z.^2;
        A(:,17) = x .* y .* (x.^2 - y.^2);
        A(:,18) = z .* y .* (3 * x.^2 - y.^2);
        A(:,19) = (6 * z.^2 - y.^2 - x.^2) .* x .* y;
        A(:,20) = (4 * z.^2 - 3 * y.^2 - 3 * x.^2) .* z .* y;
        A(:,21) = 3 * r2.^2 - 30 * r2 .* z.^2 + 35 * z.^4;
        A(:,22) = (4 * z.^2 - 3 * y.^2 - 3 * x.^2) .* x .* z;
        A(:,23) = (6 * z.^2 - y.^2 - x.^2) .* (x.^2  - y.^2);
        A(:,24) = (x.^2 - 3 * y.^2) .* x .* z;
        A(:,25) = y.^4 - 6 * x.^2 .* y.^2 + x.^4;
    end

    %----------------------------------------------------------------------
    % Order = 5 (11)
    %----------------------------------------------------------------------
    if (max_order >= 5)
        r = sqrt(x.^2 + y.^2 + z.^2);
        A(:,26) = y .* (5 * x.^4 - 10 * x.^2 .* y.^2 + y.^4);
        A(:,27) = 4 * x .* (x - y) .* y .* (x + y) .* z;
        A(:,28) = -y .* (3 * x.^2 - y.^2) .* (r - 3 * z) .* (r + 3 * z);
        A(:,29) = -2 * x .* y .* z .* (r2 - 3 * z.^2);
        A(:,30) = y .* (r2.^2 - 14 * r2 .* z.^2 + 21 * z.^4);
        A(:,31) = 1 / 2 * z .* (15 * r2.^2 - 70 * r2 .* z.^2 + 63 * z.^4);
        A(:,32) = x .* (r2.^2 - 14 * r2 .* z.^2 + 21 * z.^4);
        A(:,33) = -(x - y) .* (x + y) .* z .* (r2 - 3 * z.^2);
        A(:,34) = -x .* (x.^2 - 3 * y.^2) .* (r - 3 * z) .* (r + 3 * z);
        A(:,35) = (x.^2 - 2 * x .* y - y.^2) .* (x.^2 + 2 * x .* y - y.^2) .* z;
        A(:,36) = x .* (x.^4 - 10 * x.^2 .* y.^2 + 5 * y.^4);
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

    %% Smooth a static off-resonance map
    tic; fprintf('(%2d/%2d): Calculating a smooth approximation... ', idx3, N3);
    %----------------------------------------------------------------------
    % Calculate a smooth spherical harmonics approximation to the off-resonance map
    % Given: y = Ax
    % Ordinary least squares (OLS): ||y - Ax||_2^2
    % Weighted least squares (WLS): ||W(y - Ax)||_2^2
    %----------------------------------------------------------------------
    voxel_index = find(mask(:,:,idx3));
    im_echo_ = im_echo(:,:,idx3);
    B0map_raw_ = B0map_raw(:,:,idx3);
    w = abs(im_echo_(voxel_index));
    x_fit = bsxfun(@times, A(voxel_index,:), w) \ (w .* B0map_raw_(voxel_index));
    B0map_fit(:,:,idx3) = reshape(A * x_fit, [N1 N2]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end

%% Display results
figure('Color', 'w');
montage(reorient(B0map_fit), 'DisplayRange', []);
colormap(hot(256)); colorbar;
caxis([-60 60]);
title('Spherical harmonics fit [Hz]')
export_fig(fullfile(output_directory, 'B0map_fit'), '-r400', '-tif');

figure('Color', 'w');
montage(reorient(B0map_fit .* mask), 'DisplayRange', []);
colormap(hot(256)); colorbar;
caxis([-60 60]);
title('Masked spherical harmonics fit [Hz]')
export_fig(fullfile(output_directory, 'B0map_fit_mask'), '-r400', '-tif');

%% Fill empty space in raw B0 map with spherical harmonics fit
B0map_fill = B0map_fit;
B0map_fill(mask) = B0map_raw(mask);

figure('Color', 'w');
montage(reorient(B0map_fill), 'DisplayRange', []);
colormap(hot(256)); colorbar;
caxis([-60 60]);
title('Raw + Spherical harmonics fit [Hz]')
export_fig(fullfile(output_directory, 'B0map_fill'), '-r400', '-tif');

if 0
%% Perform TGV-based denoising on a static off-resonance map
%--------------------------------------------------------------------------
% Primal-dual method for TGV denoising
%--------------------------------------------------------------------------
B0map_tgv = zeros(N1, N2, N3, 'double');
for idx = 1:N3
    tic; fprintf('(%2d/%2d): Performing TGV smoothing (lambda = %g)... ', idx, N3, lambda);
    B0map_tgv(:,:,idx) = TGV_denoise(B0map_fill(:,:,idx), lambda, maxiter_tgv, 0);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end
end

%% Fill the holes in the mask
%----------------------------------------------------------------------
% Fill voids
%----------------------------------------------------------------------
mask_fill = bwareaopen(mask, 30);
mask_fill = imfill(mask_fill, 'holes');

%----------------------------------------------------------------------
% Dilate the mask
%----------------------------------------------------------------------
se = strel('disk', 6);
mask_fill = imdilate(mask_fill, se);

%----------------------------------------------------------------------
% Fill voids again!
%----------------------------------------------------------------------
mask_fill = imfill(mask_fill, 'holes');

%% Apply a mask on a TGV smoothed off-resonance map
%B0map_tgv_mask = B0map_tgv .* mask_fill;
B0map_fit_mask = B0map_fit .* mask_fill;

%% Display results
figure('Color', 'w');
montage(reorient(mask_fill), 'DisplayRange', []);
title('Mask without holes')
export_fig(fullfile(output_directory, 'mask_fill'), '-r400', '-tif');

if 0
figure('Color', 'w');
montage(reorient(B0map_tgv), 'DisplayRange', []);
colormap(hot(256)); colorbar;
caxis([-60 60]);
title('TGV smoothing [Hz]')
export_fig(fullfile(output_directory, 'B0map_tgv'), '-r400', '-tif');

figure('Color', 'w');
montage(reorient(B0map_tgv_mask), 'DisplayRange', []);
colormap(hot(256)); colorbar;
caxis([-60 60]);
title('TGV smoothing x mask [Hz]')
export_fig(fullfile(output_directory, 'B0map_tgv_mask'), '-r400', '-tif');
end

figure('Color', 'w');
montage(reorient(B0map_fit_mask), 'DisplayRange', []);
colormap(hot(256)); colorbar;
caxis([-60 60]);
title('Spherical harmonics fit x mask [Hz]')
export_fig(fullfile(output_directory, 'B0map_fit_mask'), '-r400', '-tif');


%     %% Display results
%     figure('Color', 'w', 'Position', [1 1 1600 823]);
%     im_montage = cat(2, flip(rot90(B0map_raw,-1),2)    , flip(rot90(B0map_raw .* voxel_mask,-1),2), ...
%                         flip(rot90(B0map_fit,-1),2)    , flip(rot90(B0map_temp,-1),2), ...
%                         flip(rot90(B0map_fit_TGV,-1),2), flip(rot90(B0map_fit_TGV .* mask,-1),2));
% %     im_montage = cat(2, B0map_raw, B0map_raw .* voxel_mask, ...
% %                         B0map_fit_TGV, B0map_fit_TGV .* voxel_mask);
%     imagesc(im_montage); axis image off;
%     caxis([-70 70]);
%     %cmap = hot(256);
%     cmap = jet(256);
%     colormap(cmap);
%     hc = colorbar;
%     title(hc, '[Hz]', 'FontSize', 14);
%     set(hc, 'FontSize', 14);
%     text(N1/2            , 0, 'Raw B0 map'             , 'Color', 'k', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
%     text(N1/2 + 1*N1, 0, 'Raw  \times voxel mask' , 'Color', 'k', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
%     text(N1/2 + 2*N1, 0, 'Spherical harmonics fit', 'Color', 'k', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
%     text(N1/2 + 3*N1, 0, 'Raw + fit'              , 'Color', 'k', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
%     text(N1/2 + 4*N1, 0, 'TGV smoothing'          , 'Color', 'k', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
%     text(N1/2 + 5*N1, 0, 'TGV  \times  mask'      , 'Color', 'k', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
%     export_fig(fullfile(user_opts_B0map.output_directory, sprintf('B0map_slice%d', idx3)), '-r400', '-tif');
%     drawnow;
% end

end