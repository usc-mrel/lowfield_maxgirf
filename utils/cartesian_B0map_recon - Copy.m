function [im_gre,B0map,x_gre,y_gre,z_gre,TE,csm] = cartesian_B0map_recon(gre_noise_fullpaths, gre_data_fullpaths, gre_dat_fullpaths, user_opts_cartesian)
% Inputs
%   gre_noise_full_paths    cell structure containing full paths to noise only h5 files
%   gre_data_full_paths     cell structure containing full paths to mutlti-echo h5 files
%   gre_dat_full_paths      cell structure containing full paths to Siemens dat files

%% Read the first ismrmrd file to get header information
start_time = tic;
tic; fprintf('Reading an ismrmrd file: %s... ', gre_data_fullpaths{1});
if exist(gre_data_fullpaths{1}, 'file')
    dset = ismrmrd.Dataset(gre_data_fullpaths{1}, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , gre_data_fullpaths{1});
end

%% Get imaging parameters from the XML header
header = ismrmrd.xml.deserialize(dset.readxml);

%--------------------------------------------------------------------------
% Encoding
%--------------------------------------------------------------------------
Nx = header.encoding.reconSpace.matrixSize.x; % number of voxels in image-space (RO,row)
Ny = header.encoding.reconSpace.matrixSize.y; % number of voxels in image-space (PE,col)
Nz = header.encoding.reconSpace.matrixSize.z; % number of voxels in image-space (SS,slice)

% Zeropad in kspace
zpad_factor = user_opts_cartesian.zpad_factor;
Nx = Nx * zpad_factor;
Ny = Ny * zpad_factor;
N  = Nx * Ny * Nz; % total number of voxels in image-space

%% Parse the ISMRMRD header
raw_data  = dset.readAcquisition(); % read all the acquisitions
nr_slices = double(max(raw_data.head.idx.slice)) + 1;

%% Set the dimension of multi-echo images
nr_TEs = length(gre_data_fullpaths);
dims = [Nx Ny max(Nz,nr_slices) nr_TEs];

%% Preallocate memory
im_gre = complex(zeros(dims, 'single'));
TE     = zeros(nr_TEs,1, 'single'); % [sec]

%% Reconstruct multi-echo GRE data
for idx = 1:nr_TEs
    %----------------------------------------------------------------------
    % Perform Cartesian reconstruction
    %----------------------------------------------------------------------
    if idx == 1
        [im,header,x,y,z,csm] = cartesian_recon(gre_noise_fullpaths{idx}, gre_data_fullpaths{idx}, gre_dat_fullpaths{idx}, user_opts_cartesian);
        csm = squeeze(csm);
    else
        [im,header,x,y,z] = cartesian_recon(gre_noise_fullpaths{idx}, gre_data_fullpaths{idx}, gre_dat_fullpaths{idx}, user_opts_cartesian);
    end

    %----------------------------------------------------------------------
    % Save the reconstructed image
    %----------------------------------------------------------------------
    im_gre(:,:,:,idx) = reshape(mean(im,4), dims(1:3));
    TE(idx) = header.sequenceParameters.TE * 1e-3; % [msec] * [sec/1e3msec] => [sec]
end

%% Get spatial coordinates in DCS [mm]
x_gre = reshape(x, dims(1:3)); % Nx x Ny x Nz
y_gre = reshape(y, dims(1:3)); % Nx x Ny x Nz
z_gre = reshape(z, dims(1:3)); % Nx x Ny x Nz

%% Calculate a smoothed off-resonance map per slice
start_time = tic;
B0map = zeros(dims(1:3), 'single');

for idx3 = 1:dims(3)
    %% Calculate a raw static off-resonance map [Hz]
    tic; fprintf('(%2d/%2d): Calculating a static off-resonance map... ', idx3, dims(3));
    B0map_raw = zeros(dims(1:2), 'single'); % [Hz]
    nr_voxels = prod(dims(1:2));
    for idx = 1:nr_voxels
        %------------------------------------------------------------------
        % Perform least-squares fitting per voxel
        % m = a * TE + b, m in [rad], a in [rad/sec], TE in [sec]
        % m(1) = a * TE(1) + b    [m(1)] = [TE(1) 1][a]
        % m(2) = a * TE(2) + b => [m(2)] = [TE(2) 1][b] => m = Av
        % m(N) = a * TE(N) + b    [m(N)] = [TE(N) 1]
        %------------------------------------------------------------------
        [idx1,idx2] = ind2sub(dims(1:2), idx);
        m = unwrap(reshape(angle(im_gre(idx1,idx2,idx3,:)), [dims(4) 1])); % [rad]
        A = [TE ones(dims(4),1)]; % [sec]
        v = A \ m; % [rad/sec]
        B0map_raw(idx) = v(1) / (2 * pi); % [rad/sec] * [cycle/2pi rad] => [Hz]
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

    %% Set spatial coordinates of the current slice
    x = reshape(x_gre(:,:,idx3), [N 1]); % Nx x Ny => N x 1
    y = reshape(y_gre(:,:,idx3), [N 1]); % Nx x Ny => N x 1
    z = reshape(z_gre(:,:,idx3), [N 1]); % Nx x Ny => N x 1

    %% Calculate real-valued spherical harmonic basis functions
    tic; fprintf('(%2d/%2d): Calculating real-valued spherical harmonic basis functions... ', idx3, dims(3));
    %----------------------------------------------------------------------
    % Calculate the number of basis functions
    %----------------------------------------------------------------------
    max_order = 5;
    Ns = 0; % number of basis functions
    for order = 0:max_order
        Ns = Ns + 2 * order + 1;
    end
    A = zeros(N, Ns, 'single');

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

    %% Calculate a voxel mask
    im_ = im_gre(:,:,idx3,1); % first TE image
    % 0.4 for a phantom
    voxel_mask = (abs(im_) > 0.4 * mean(abs(reshape(im_, [N 1]))));
    voxel_mask = bwareaopen(voxel_mask, 2);

    %----------------------------------------------------------------------
    % Fill voids
    %----------------------------------------------------------------------
    mask = bwareaopen(voxel_mask, 30);
    mask = imfill(mask, 'holes');

    %----------------------------------------------------------------------
    % Dilate the mask
    %----------------------------------------------------------------------
    se = strel('disk', 5);
    mask = imdilate(mask, se);

    %----------------------------------------------------------------------
    % Fill voids again!
    %----------------------------------------------------------------------
    mask = imfill(mask, 'holes');
    voxel_mask = (voxel_mask .* mask > 0);

    %% Smooth a static off-resonance map
    tic; fprintf('(%2d/%2d): Calculating a smooth approximation... ', idx3, dims(3));
    %----------------------------------------------------------------------
    % Calculate a smooth spherical harmonics approximation to the off-resonance map
    % Given: y = Ax
    % Ordinary least squares (OLS): ||y - Ax||_2^2
    % Weighted least squares (WLS): ||W(y - Ax)||_2^2
    %----------------------------------------------------------------------
    voxel_index = find(voxel_mask);
    w = abs(im_(voxel_index));
    x_fit = bsxfun(@times, A(voxel_index,:), w) \ (w .* B0map_raw(voxel_index));
    B0map_fit = reshape(A * x_fit, dims(1:2));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

    %----------------------------------------------------------------------
    % Primal-dual method for TGV denoising
    %----------------------------------------------------------------------
    tic; fprintf('(%2d/%2d): Performing TGV smoothing... ', idx3, dims(3));
    %B0map_temp = B0map_fit;
    %B0map_temp(voxel_mask) = B0map_raw_(voxel_mask);
    B0map_temp = B0map_raw .* voxel_mask;
    maxiter_tgv = 1500;
    lambda = 1e2; % regularization parameter % human?
    %lambda = 1e2; % regularization parameter
    B0map_fit_TGV = TGV_denoise(B0map_temp, lambda, maxiter_tgv, 0);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

    %% Get a smoothed static off-resonance map
    B0map(:,:,idx3) = B0map_fit_TGV .* mask;

    %% Display results
    figure('Color', 'w', 'Position', [1 1 1600 823]);
    im_montage = cat(2, flip(rot90(B0map_raw,-1),2)    , flip(rot90(B0map_raw .* voxel_mask,-1),2), ...
                        flip(rot90(B0map_fit_TGV,-1),2), flip(rot90(B0map_fit_TGV .* voxel_mask,-1),2));
%     im_montage = cat(2, B0map_raw, B0map_raw .* voxel_mask, ...
%                         B0map_fit_TGV, B0map_fit_TGV .* voxel_mask);
    imagesc(im_montage); axis image off;
    caxis([-70 70]);
    cmap = jet(256);
    colormap(cmap);
    hc = colorbar;
    title(hc, '[Hz]', 'FontSize', 14);
    set(hc, 'FontSize', 14);
    text(dims(1)/2            , 0, 'Raw B0 map'                     , 'Color', 'k', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    text(dims(1)/2 + 1*dims(1), 0, 'Raw B0 map \times voxel mask'   , 'Color', 'k', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    text(dims(1)/2 + 2*dims(1), 0, 'TGV smoothing'                  , 'Color', 'k', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    text(dims(1)/2 + 3*dims(1), 0, 'TGV smoothing \times voxel mask', 'Color', 'k', 'FontSize', 14, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    export_fig(sprintf('B0map_slice%d', idx3), '-r400', '-tif');
    drawnow;
end

end