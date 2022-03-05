% demo_cartesian_recon_GRE_datasets.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/16/2022, Last modified: 01/18/2022

%% Clean slate
close all; clear all; clc;

%% Set source directories
src_directory = 'E:\lowfield_maxgirf';
ismrmrd_directory = 'D:\ismrmrd\ismrmrd';

%% Add source directories to search path
addpath(genpath(src_directory));
addpath(genpath(ismrmrd_directory));

%% Define data directory
%--------------------------------------------------------------------------
% noise only ISMRMRD format
%--------------------------------------------------------------------------
ismrmrd_noise_fullpath{1} = 'D:\lowfield\NHLBI\data\20201102_NV_brain\noise\noise_meas_MID00261_FID03644_gre_TE1.h5';
ismrmrd_noise_fullpath{2} = 'D:\lowfield\NHLBI\data\20201102_NV_brain\noise\noise_meas_MID00262_FID03645_gre_TE2.h5';
ismrmrd_noise_fullpath{3} = 'D:\lowfield\NHLBI\data\20201102_NV_brain\noise\noise_meas_MID00263_FID03646_gre_TE3.h5';
ismrmrd_noise_fullpath{4} = 'D:\lowfield\NHLBI\data\20201102_NV_brain\noise\noise_meas_MID00264_FID03647_gre_TE4.h5';
ismrmrd_noise_fullpath{5} = 'D:\lowfield\NHLBI\data\20201102_NV_brain\noise\noise_meas_MID00265_FID03648_gre_TE5.h5';

%--------------------------------------------------------------------------
% ISMRMRD format
%--------------------------------------------------------------------------
ismrmrd_data_fullpath{1}  = 'D:\lowfield\NHLBI\data\20201102_NV_brain\h5\meas_MID00261_FID03644_gre_TE1.h5';
ismrmrd_data_fullpath{2}  = 'D:\lowfield\NHLBI\data\20201102_NV_brain\h5\meas_MID00262_FID03645_gre_TE2.h5';
ismrmrd_data_fullpath{3}  = 'D:\lowfield\NHLBI\data\20201102_NV_brain\h5\meas_MID00263_FID03646_gre_TE3.h5';
ismrmrd_data_fullpath{4}  = 'D:\lowfield\NHLBI\data\20201102_NV_brain\h5\meas_MID00264_FID03647_gre_TE4.h5';
ismrmrd_data_fullpath{5}  = 'D:\lowfield\NHLBI\data\20201102_NV_brain\h5\meas_MID00265_FID03648_gre_TE5.h5';

%--------------------------------------------------------------------------
% Siemens raw data format
%--------------------------------------------------------------------------
siemens_dat_fullpath{1}  = 'D:\lowfield\NHLBI\data\20201102_NV_brain\meas_MID00261_FID03644_gre_TE1.dat';
siemens_dat_fullpath{2}  = 'D:\lowfield\NHLBI\data\20201102_NV_brain\meas_MID00262_FID03645_gre_TE2.dat';
siemens_dat_fullpath{3}  = 'D:\lowfield\NHLBI\data\20201102_NV_brain\meas_MID00263_FID03646_gre_TE3.dat';
siemens_dat_fullpath{4}  = 'D:\lowfield\NHLBI\data\20201102_NV_brain\meas_MID00264_FID03647_gre_TE4.dat';
siemens_dat_fullpath{5}  = 'D:\lowfield\NHLBI\data\20201102_NV_brain\meas_MID00265_FID03648_gre_TE5.dat';
Ne = length(ismrmrd_noise_fullpath);

% "phase_sign" and "read_sign" can be determined only from Siemens raw data 
% format now until the ISMRMRD format includes these as part of its header
user_opts.phase_sign = 1;
user_opts.read_sign = -1;
user_opts.remove_oversampling = 0; % Remove readout oversampling

save_filename = sprintf('human_sagittal_data_ro%d', user_opts.remove_oversampling);

%% Read the first echo dataset to calculate coil sensitivity maps
[kspace,header,r_dcs] = siemens_read_kspace_cartesian_multislice(ismrmrd_noise_fullpath{1}, ismrmrd_data_fullpath{1}, [], user_opts);
[N1,N2,Ns,Nc] = size(kspace);

%% Zeropad in k-space
N1_zpad = 2 * N1;
N2_zpad = 2 * N2;
kspace = zpad(kspace, [N1_zpad N2_zpad Ns Nc]);

%% Calculate coil sensitivity maps
start_time = tic;
csm = complex(zeros(N1_zpad, N2_zpad, Ns, Nc, 'double'));
for idx = 1:Ns
    tstart = tic; fprintf('(%2d/%2d): Calculating coil sensitivity maps with Walsh method... ', idx, Ns);
    %----------------------------------------------------------------------
    % Calculate the calibration region of k-space
    %----------------------------------------------------------------------
    cal_shape = [32 32];
    cal_data = crop(reshape(kspace(:,:,idx,:), [N1_zpad N2_zpad Nc]), [cal_shape Nc]);
    cal_data = bsxfun(@times, cal_data, hamming(cal_shape(1)) * hamming(cal_shape(2)).');

    %----------------------------------------------------------------------
    % Calculate coil sensitivity maps
    %----------------------------------------------------------------------
    cal_im = zpad(cal_data, [N1_zpad N2_zpad Nc]);
    for dim = 1:2
        cal_im = 1 / sqrt(size(cal_im,dim)) * fftshift(fft(ifftshift(cal_im, dim), [], dim), dim);
    end
    csm(:,:,idx,:) = ismrm_estimate_csm_walsh(cal_im);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Perform FFT reconstruction for all echo datasets
im_echo = complex(zeros(N1_zpad, N2_zpad, Ns, Ne, 'double'));
TEs = zeros(Ne, 1, 'double'); % echo time in [sec]
for idx = 1:Ne
    %----------------------------------------------------------------------
    % Read in k-space data
    %----------------------------------------------------------------------
    [kspace,header] = siemens_read_kspace_cartesian_multislice(ismrmrd_noise_fullpath{idx}, ismrmrd_data_fullpath{idx}, [], user_opts); % N1 x N2 x Ns x Nc

    %----------------------------------------------------------------------
    % Zeropad in k-space
    %----------------------------------------------------------------------
    kspace = zpad(kspace, [N1_zpad N2_zpad Ns Nc]);

    %----------------------------------------------------------------------
    % Perform FFT (k-space <=> image-space)
    %----------------------------------------------------------------------
    imc_fft = kspace;
    for dim = 1:2
        imc_fft = 1 / sqrt(size(imc_fft,dim)) * fftshift(fft(ifftshift(imc_fft, dim), [], dim), dim); % N1_zpad x N2_zpad x Ns x Nc
    end

    %----------------------------------------------------------------------
    % Perform optimal coil combination 
    %----------------------------------------------------------------------
    im_echo(:,:,:,idx) = sum(bsxfun(@times, conj(csm), imc_fft), 4); % N1_zpad x N2_zpad x Ns

    %----------------------------------------------------------------------
    % Collect TE
    %----------------------------------------------------------------------
    TEs(idx) = header.sequenceParameters.TE * 1e-3; % [msec] * [sec/1e3msec] => [sec]
end

%% Save results
save(save_filename, 'im_echo', 'TEs', '-v7.3');

%% Display images
for idx = 1:Ns
    slice_nr = idx;
    %----------------------------------------------------------------------
    % Display the magnitude of GRE images
    %----------------------------------------------------------------------
    figure('Color', 'k', 'Position', [4 201 1553 597]);
    color_order = get(gca, 'colororder');
    imagesc(reshape(abs(im_echo(:,:,slice_nr,:)), [N1_zpad N2_zpad*Ne])); axis image;
    set(gca, 'XColor', 'w', 'YColor', 'w', 'XTick', [], 'YTick', []);
    caxis([0 16/2]);
    colormap(gray(256));
    hc = colorbar;
    set(hc, 'Color', 'w', 'FontSize', 14);
    for idx1 = 1:Ne
        text(N2_zpad/2 + (idx1 - 1) * N2_zpad, 0, sprintf('TE = %4.2f ms', TEs(idx1)*1e3), 'Color', 'w', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
    end
    text(N2_zpad / 2 + N2_zpad * floor(Ne/2), 0, 'Magnitude of Cartesian GRE images', 'FontSize', 20, 'Color', color_order(3,:), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    export_fig(sprintf('GRE_images_slice%d_mag', slice_nr), '-r300', '-tif', '-c[70,380,200,560]'); % [top,right,bottom,left]
    close gcf;
 
    %----------------------------------------------------------------------
    % Display the phase of GRE images
    %----------------------------------------------------------------------
    figure('Color', 'k', 'Position', [4 201 1553 597]);
    imagesc(reshape(angle(im_echo(:,:,slice_nr,:))*180/pi, [N1_zpad N2_zpad*Ne])); axis image;
    set(gca, 'XColor', 'w', 'YColor', 'w', 'XTick', [], 'YTick', []);
    colormap(hsv(256));
    hc = colorbar;
    set(hc, 'Color', 'w', 'FontSize', 14);
    for idx1 = 1:Ne
        text(N2_zpad/2 + (idx1 - 1) * N2_zpad, 0, sprintf('TE = %4.2f ms', TEs(idx1)*1e3), 'Color', 'k', 'FontSize', 14, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
    end
    text(N2_zpad / 2 + N2_zpad * floor(Ne/2), 0, 'Phase of Cartesian GRE images', 'FontSize', 20, 'Color', color_order(3,:), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
    export_fig(sprintf('GRE_images_slice%d_phase', slice_nr), '-r300', '-tif', '-c[70,380,200,560]'); % [top,right,bottom,left]
    close gcf;
end