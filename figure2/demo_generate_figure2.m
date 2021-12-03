% demo_generate_figure2.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 05/20/2021, Last modified: 05/20/2021

%% Clean slate
close all; clear; clc;

%% Add paths
thirdparty_directory = 'D:\lowfield_maxgirf\thirdparty';
addpath(genpath(thirdparty_directory));

%% Set the directory names
axial_data_directory1 = 'D:\lowfield_maxgirf\figure2\20201102_NV_brain_se_spiral_1102_ax_s24_osf1_B0_correction1_concomitant_correction0_Lmax50_L8';
axial_data_directory2 = 'D:\lowfield_maxgirf\figure2\20201102_NV_brain_se_spiral_1102_ax_s24_osf1_B0_correction0_concomitant_correction1_Lmax50_L8';
axial_data_directory3 = 'D:\lowfield_maxgirf\figure2\20201102_NV_brain_se_spiral_1102_ax_s24_osf1_B0_correction1_concomitant_correction1_Lmax50_L8';

sagittal_data_directory1 = 'D:\lowfield_maxgirf\figure2\20201102_NV_brain_se_spiral_1102_sag_s24_osf2_B0_correction1_concomitant_correction0_Lmax100_L30';
sagittal_data_directory2 = 'D:\lowfield_maxgirf\figure2\20201102_NV_brain_se_spiral_1102_sag_s24_osf2_B0_correction0_concomitant_correction1_Lmax100_L30';
sagittal_data_directory3 = 'D:\lowfield_maxgirf\figure2\20201102_NV_brain_se_spiral_1102_sag_s24_osf2_B0_correction1_concomitant_correction1_Lmax100_L30';

%% Calculate the NRMSE between full-rank and rank-L reconstructions
slice_offsets_axial = [-70; -35; 0; 35; 70; 105; -52.5; -17.5; 17.5; 52.5; 87.5] * 1e-3; % [m]
actual_slice_nrs_axial = [1,3,5,7,9,11,2,4,6,8,10].';
Lmax_axial = 50;
N3 = 11;

NRMSE_axial1 = zeros(Lmax_axial, N3, 'double');
NRMSE_axial2 = zeros(Lmax_axial, N3, 'double');
NRMSE_axial3 = zeros(Lmax_axial, N3, 'double');

start_time = tic;
for idx2 = 1:N3
    tstart = tic; fprintf('Processing slice (%d/%d)... ', idx2, N3);
    %----------------------------------------------------------------------
    % Load mat files
    %----------------------------------------------------------------------
    maxgirf_cpr_filename = sprintf('maxgirf_cpr_slice%d.mat', idx2);
    load(fullfile(axial_data_directory1, maxgirf_cpr_filename)); % B0=1, conc=0
    im_maxgirf_cpr1 = im_maxgirf_cpr;
    load(fullfile(axial_data_directory2, maxgirf_cpr_filename)); % B0=0, conc=1
    im_maxgirf_cpr2 = im_maxgirf_cpr;
    load(fullfile(axial_data_directory3, maxgirf_cpr_filename)); % B0=1, conc=1
    im_maxgirf_cpr3 = im_maxgirf_cpr;

    %----------------------------------------------------------------------
    % Calculate a mask defining the image and noise regions using 
    % the iterative intermean algorithm
    %----------------------------------------------------------------------
    m_full_axial1 = im_maxgirf_cpr1(:,:,Lmax_axial);
    m_full_axial2 = im_maxgirf_cpr2(:,:,Lmax_axial);
    m_full_axial3 = im_maxgirf_cpr3(:,:,Lmax_axial);
    [N1,N2] = size(m_full_axial3);

    %----------------------------------------------------------------------
    % Calculate the maximum magnitude
    %----------------------------------------------------------------------
    mask_espirit = abs(m_full_axial3) > 0;
    im_max = abs(m_full_axial3(mask_espirit));
    figure, imagesc(mask_espirit); axis image; drawnow;

    %----------------------------------------------------------------------
    % Calculate a mask
    %----------------------------------------------------------------------
    level = isodata(im_max, 'log');
    mask = false(N1, N2);
    mask(abs(m_full_axial3) > level) = 1;
    figure, imagesc(mask); axis image; drawnow;

    %----------------------------------------------------------------------
    % Calculate the NRMSE between full-rank and rank-L reconstructions
    %----------------------------------------------------------------------
    for idx1 = 1:Lmax_axial
        m_L = im_maxgirf_cpr1(:,:,idx1);
        NRMSE_axial1(idx1,idx2) = norm((m_full_axial1(:) - m_L(:)) .* mask(:)) / norm(m_full_axial1(:) .* mask(:));

        m_L = im_maxgirf_cpr2(:,:,idx1);
        NRMSE_axial2(idx1,idx2) = norm((m_full_axial2(:) - m_L(:)) .* mask(:)) / norm(m_full_axial2(:) .* mask(:));

        m_L = im_maxgirf_cpr3(:,:,idx1);
        NRMSE_axial3(idx1,idx2) = norm((m_full_axial3(:) - m_L(:)) .* mask(:)) / norm(m_full_axial3(:) .* mask(:));
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Calculate the NRMSE between full-rank and rank-L reconstructions
slice_offsets_sagittal = [-62.47249;
                          -37.48349;
                          -12.49450;
                           12.49450;
                           37.48349;
                           62.47249;
                          -49.97799;
                          -24.98900;
                           0;
                           24.98900;
                           49.97799] * 1e-3; % [m]
actual_slice_nrs_sagittal = [1,3,5,7,9,11,2,4,6,8,10].';
Lmax_sagittal = 100;
N3 = 11;

NRMSE_sagittal1 = zeros(Lmax_sagittal, N3, 'double');
NRMSE_sagittal2 = zeros(Lmax_sagittal, N3, 'double');
NRMSE_sagittal3 = zeros(Lmax_sagittal, N3, 'double');

start_time = tic;
for idx2 = 1:N3
    tstart = tic; fprintf('Processing slice (%d/%d)... ', idx2, N3);
    %----------------------------------------------------------------------
    % Load mat files
    %----------------------------------------------------------------------
    maxgirf_cpr_filename = sprintf('maxgirf_cpr_slice%d.mat', idx2);
    load(fullfile(sagittal_data_directory1, maxgirf_cpr_filename)); % B0=1, conc=0
    im_maxgirf_cpr1 = im_maxgirf_cpr;
    load(fullfile(sagittal_data_directory2, maxgirf_cpr_filename)); % B0=0, conc=1
    im_maxgirf_cpr2 = im_maxgirf_cpr;
    load(fullfile(sagittal_data_directory3, maxgirf_cpr_filename)); % B0=1, conc=1
    im_maxgirf_cpr3 = im_maxgirf_cpr;

    %----------------------------------------------------------------------
    % Calculate a mask defining the image and noise regions using 
    % the iterative intermean algorithm
    %----------------------------------------------------------------------
    m_full_sagittal1 = im_maxgirf_cpr1(:,:,Lmax_sagittal);
    m_full_sagittal2 = im_maxgirf_cpr2(:,:,Lmax_sagittal);
    m_full_sagittal3 = im_maxgirf_cpr3(:,:,Lmax_sagittal);
    [N1,N2] = size(m_full_sagittal3);

    %----------------------------------------------------------------------
    % Calculate the maximum magnitude
    %----------------------------------------------------------------------
    mask_espirit = abs(m_full_sagittal3) > 0;
    im_max = abs(m_full_sagittal3(mask_espirit));
    figure, imagesc(mask_espirit); axis image; drawnow;

    %----------------------------------------------------------------------
    % Calculate a mask
    %----------------------------------------------------------------------
    level = isodata(im_max, 'log');
    mask = false(N1, N2);
    mask(abs(m_full_sagittal3) > level) = 1;
    figure, imagesc(mask); axis image; drawnow;

    %----------------------------------------------------------------------
    % Calculate the NRMSE between full-rank and rank-L reconstructions
    %----------------------------------------------------------------------
    for idx1 = 1:Lmax_sagittal
        m_L = im_maxgirf_cpr1(:,:,idx1);
        NRMSE_sagittal1(idx1,idx2) = norm((m_full_sagittal1(:) - m_L(:)) .* mask(:)) / norm(m_full_sagittal1(:) .* mask(:));

        m_L = im_maxgirf_cpr2(:,:,idx1);
        NRMSE_sagittal2(idx1,idx2) = norm((m_full_sagittal2(:) - m_L(:)) .* mask(:)) / norm(m_full_sagittal2(:) .* mask(:));

        m_L = im_maxgirf_cpr3(:,:,idx1);
        NRMSE_sagittal3(idx1,idx2) = norm((m_full_sagittal3(:) - m_L(:)) .* mask(:)) / norm(m_full_sagittal3(:) .* mask(:));
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Prepare inset axial images 
%--------------------------------------------------------------------------
% Load a .mat file
%--------------------------------------------------------------------------
maxgirf_cpr_filename = sprintf('maxgirf_cpr_slice%d.mat', 9);
load(fullfile(axial_data_directory1, maxgirf_cpr_filename)); % B0=1, conc=0
im_maxgirf_cpr1 = im_maxgirf_cpr;
load(fullfile(axial_data_directory2, maxgirf_cpr_filename)); % B0=0, conc=1
im_maxgirf_cpr2 = im_maxgirf_cpr;
load(fullfile(axial_data_directory3, maxgirf_cpr_filename)); % B0=1, conc=1
im_maxgirf_cpr3 = im_maxgirf_cpr;

%--------------------------------------------------------------------------
% Calculate a mask defining the image and noise regions using 
% the iterative intermean algorithm
%--------------------------------------------------------------------------
m_full_axial1 = im_maxgirf_cpr1(:,:,Lmax_axial);
m_full_axial2 = im_maxgirf_cpr2(:,:,Lmax_axial);
m_full_axial3 = im_maxgirf_cpr3(:,:,Lmax_axial);
[N1,N2] = size(m_full_axial3);

%--------------------------------------------------------------------------
% Calculate the maximum magnitude
%--------------------------------------------------------------------------
mask_espirit = abs(m_full_axial3) > 0;
im_max = abs(m_full_axial3(mask_espirit));

%--------------------------------------------------------------------------
% Calculate a mask
%--------------------------------------------------------------------------
level = isodata(im_max, 'log');
mask = false(N1, N2);
mask(abs(m_full_axial3) > level) = 1;

%--------------------------------------------------------------------------
% Inset image
%--------------------------------------------------------------------------
reorient = @(x) flip(rot90(x, -1), 2);

L1_axial = 3;
L2_axial = 4;
L3_axial = 6;
L4_axial = 8;

idx1_range = (40:280).';
idx2_range = (30:320).';
N1_zoom1 = length(idx1_range);
N2_zoom1 = length(idx2_range);
idx2_range = (-floor(271/2):ceil(271/2)-1).' + floor(N1/2) + 1  + 15;
idx1_range = (-floor(251/2):ceil(251/2)-1).' + floor(N2/2) + 1;
N1_zoom1 = length(idx1_range);
N2_zoom1 = length(idx2_range);

im_diff1_1 = (m_full_axial1 - im_maxgirf_cpr1(:,:,L1_axial)) ./ abs(m_full_axial1) * 100 .* mask;
im_diff1_2 = (m_full_axial1 - im_maxgirf_cpr1(:,:,L2_axial)) ./ abs(m_full_axial1) * 100 .* mask;
im_diff1_3 = (m_full_axial1 - im_maxgirf_cpr1(:,:,L3_axial)) ./ abs(m_full_axial1) * 100 .* mask;
im_diff1_4 = (m_full_axial1 - im_maxgirf_cpr1(:,:,L4_axial)) ./ abs(m_full_axial1) * 100 .* mask;

im_diff2_1 = (m_full_axial2 - im_maxgirf_cpr2(:,:,L1_axial)) ./ abs(m_full_axial2) * 100 .* mask;
im_diff2_2 = (m_full_axial2 - im_maxgirf_cpr2(:,:,L2_axial)) ./ abs(m_full_axial2) * 100 .* mask;
im_diff2_3 = (m_full_axial2 - im_maxgirf_cpr2(:,:,L3_axial)) ./ abs(m_full_axial2) * 100 .* mask;
im_diff2_4 = (m_full_axial2 - im_maxgirf_cpr2(:,:,L4_axial)) ./ abs(m_full_axial2) * 100 .* mask;

im_diff3_1 = (m_full_axial3 - im_maxgirf_cpr3(:,:,L1_axial)) ./ abs(m_full_axial3) * 100 .* mask;
im_diff3_2 = (m_full_axial3 - im_maxgirf_cpr3(:,:,L2_axial)) ./ abs(m_full_axial3) * 100 .* mask;
im_diff3_3 = (m_full_axial3 - im_maxgirf_cpr3(:,:,L3_axial)) ./ abs(m_full_axial3) * 100 .* mask;
im_diff3_4 = (m_full_axial3 - im_maxgirf_cpr3(:,:,L4_axial)) ./ abs(m_full_axial3) * 100 .* mask;

nan_mask = zeros(N1,N2);
nan_mask(mask) = 1;
nan_mask(~mask) = NaN;

im_diff1_1 = im_diff1_1 .* nan_mask;
im_diff1_2 = im_diff1_2 .* nan_mask;
im_diff1_3 = im_diff1_3 .* nan_mask;
im_diff1_4 = im_diff1_4 .* nan_mask;

im_diff2_1 = im_diff2_1 .* nan_mask;
im_diff2_2 = im_diff2_2 .* nan_mask;
im_diff2_3 = im_diff2_3 .* nan_mask;
im_diff2_4 = im_diff2_4 .* nan_mask;

im_diff3_1 = im_diff3_1 .* nan_mask;
im_diff3_2 = im_diff3_2 .* nan_mask;
im_diff3_3 = im_diff3_3 .* nan_mask;
im_diff3_4 = im_diff3_4 .* nan_mask;

diff_montage_axial1 = cat(1, cat(2, reorient(im_diff1_1(idx1_range, idx2_range)), reorient(im_diff1_2(idx1_range, idx2_range))), ...
                             cat(2, reorient(im_diff1_3(idx1_range, idx2_range)), reorient(im_diff1_4(idx1_range, idx2_range))));

diff_montage_axial2 = cat(1, cat(2, reorient(im_diff2_1(idx1_range, idx2_range)), reorient(im_diff2_2(idx1_range, idx2_range))), ...
                             cat(2, reorient(im_diff2_3(idx1_range, idx2_range)), reorient(im_diff2_4(idx1_range, idx2_range))));

diff_montage_axial3 = cat(1, cat(2, reorient(im_diff3_1(idx1_range, idx2_range)), reorient(im_diff3_2(idx1_range, idx2_range))), ...
                             cat(2, reorient(im_diff3_3(idx1_range, idx2_range)), reorient(im_diff3_4(idx1_range, idx2_range))));

%% Prepare inset sagittal images 
%--------------------------------------------------------------------------
% Load a .mat file
%--------------------------------------------------------------------------
maxgirf_cpr_filename = sprintf('maxgirf_cpr_slice%d.mat', 11);
load(fullfile(sagittal_data_directory1, maxgirf_cpr_filename)); % B0=1, conc=0
im_maxgirf_cpr1 = im_maxgirf_cpr;
load(fullfile(sagittal_data_directory2, maxgirf_cpr_filename)); % B0=0, conc=1
im_maxgirf_cpr2 = im_maxgirf_cpr;
load(fullfile(sagittal_data_directory3, maxgirf_cpr_filename)); % B0=1, conc=1
im_maxgirf_cpr3 = im_maxgirf_cpr;

%--------------------------------------------------------------------------
% Calculate a mask defining the image and noise regions using 
% the iterative intermean algorithm
%--------------------------------------------------------------------------
m_full_sagittal1 = im_maxgirf_cpr1(:,:,Lmax_sagittal);
m_full_sagittal2 = im_maxgirf_cpr2(:,:,Lmax_sagittal);
m_full_sagittal3 = im_maxgirf_cpr3(:,:,Lmax_sagittal);
[N1,N2] = size(m_full_sagittal3);

%--------------------------------------------------------------------------
% Calculate the maximum magnitude
%--------------------------------------------------------------------------
mask_espirit = abs(m_full_sagittal3) > 0;
im_max = abs(m_full_sagittal3(mask_espirit));

%--------------------------------------------------------------------------
% Calculate a mask
%--------------------------------------------------------------------------
level = isodata(im_max, 'log');
mask = false(N1, N2);
mask(abs(m_full_sagittal3) > level) = 1;

%--------------------------------------------------------------------------
% Inset image
%--------------------------------------------------------------------------
reorient = @(x) x;

idx1_range = (180:450).';
idx2_range = (200:450).';
N1_zoom2 = length(idx1_range);
N2_zoom2 = length(idx2_range);

im_diff1_1 = (m_full_sagittal1 - im_maxgirf_cpr1(:,:,15)) ./ abs(m_full_sagittal1) * 100 .* mask;
im_diff1_2 = (m_full_sagittal1 - im_maxgirf_cpr1(:,:,20)) ./ abs(m_full_sagittal1) * 100 .* mask;
im_diff1_3 = (m_full_sagittal1 - im_maxgirf_cpr1(:,:,30)) ./ abs(m_full_sagittal1) * 100 .* mask;
im_diff1_4 = (m_full_sagittal1 - im_maxgirf_cpr1(:,:,50)) ./ abs(m_full_sagittal1) * 100 .* mask;

im_diff2_1 = (m_full_sagittal2 - im_maxgirf_cpr2(:,:,15)) ./ abs(m_full_sagittal2) * 100 .* mask;
im_diff2_2 = (m_full_sagittal2 - im_maxgirf_cpr2(:,:,20)) ./ abs(m_full_sagittal2) * 100 .* mask;
im_diff2_3 = (m_full_sagittal2 - im_maxgirf_cpr2(:,:,30)) ./ abs(m_full_sagittal2) * 100 .* mask;
im_diff2_4 = (m_full_sagittal2 - im_maxgirf_cpr2(:,:,50)) ./ abs(m_full_sagittal2) * 100 .* mask;

im_diff3_1 = (m_full_sagittal3 - im_maxgirf_cpr3(:,:,15)) ./ abs(m_full_sagittal3) * 100 .* mask;
im_diff3_2 = (m_full_sagittal3 - im_maxgirf_cpr3(:,:,20)) ./ abs(m_full_sagittal3) * 100 .* mask;
im_diff3_3 = (m_full_sagittal3 - im_maxgirf_cpr3(:,:,30)) ./ abs(m_full_sagittal3) * 100 .* mask;
im_diff3_4 = (m_full_sagittal3 - im_maxgirf_cpr3(:,:,50)) ./ abs(m_full_sagittal3) * 100 .* mask;

nan_mask = zeros(N1,N2);
nan_mask(mask) = 1;
nan_mask(~mask) = NaN;

im_diff1_1 = im_diff1_1 .* nan_mask;
im_diff1_2 = im_diff1_2 .* nan_mask;
im_diff1_3 = im_diff1_3 .* nan_mask;
im_diff1_4 = im_diff1_4 .* nan_mask;

im_diff2_1 = im_diff2_1 .* nan_mask;
im_diff2_2 = im_diff2_2 .* nan_mask;
im_diff2_3 = im_diff2_3 .* nan_mask;
im_diff2_4 = im_diff2_4 .* nan_mask;

im_diff3_1 = im_diff3_1 .* nan_mask;
im_diff3_2 = im_diff3_2 .* nan_mask;
im_diff3_3 = im_diff3_3 .* nan_mask;
im_diff3_4 = im_diff3_4 .* nan_mask;

diff_montage_sagittal1 = cat(1, cat(2, reorient(im_diff1_1(idx1_range, idx2_range)), reorient(im_diff1_2(idx1_range, idx2_range))), ...
                                cat(2, reorient(im_diff1_3(idx1_range, idx2_range)), reorient(im_diff1_4(idx1_range, idx2_range))));

diff_montage_sagittal2 = cat(1, cat(2, reorient(im_diff2_1(idx1_range, idx2_range)), reorient(im_diff2_2(idx1_range, idx2_range))), ...
                                cat(2, reorient(im_diff2_3(idx1_range, idx2_range)), reorient(im_diff2_4(idx1_range, idx2_range))));

diff_montage_sagittal3 = cat(1, cat(2, reorient(im_diff3_1(idx1_range, idx2_range)), reorient(im_diff3_2(idx1_range, idx2_range))), ...
                                cat(2, reorient(im_diff3_3(idx1_range, idx2_range)), reorient(im_diff3_4(idx1_range, idx2_range))));

%% Display Figure 2
[~,sorted_index_axial] = sort(slice_offsets_axial);
[~,sorted_index_sagittal] = sort(slice_offsets_sagittal);

%figure('Color', 'w', 'Position', [3 213 849 600]);
figure('Color', 'w', 'Position', [3 213 1032 600]);
color_order = get(gca, 'colororder');

%--------------------------------------------------------------------------
% Axial, B0=1, conc=0
%--------------------------------------------------------------------------
ax1 = subplot(2,3,1); grid on; grid minor;
hold on;
cmap = parula(7);
idx_list = [1; 2; 3; 9; 10; 11];
for idx = 1:length(idx_list)
    plot(1:Lmax_axial, NRMSE_axial1(:,idx_list(idx))*100, 'LineWidth', 1, 'Color', cmap(7-idx,:));
end
plot(L1_axial, NRMSE_axial1(L1_axial,9)*100, '.', 'LineWidth', 1, 'Color', cmap(7-4,:), 'MarkerSize', 15);
plot(L2_axial, NRMSE_axial1(L2_axial,9)*100, '.', 'LineWidth', 1, 'Color', cmap(7-4,:), 'MarkerSize', 15);
plot(L3_axial, NRMSE_axial1(L3_axial,9)*100, '.', 'LineWidth', 1, 'Color', cmap(7-4,:), 'MarkerSize', 15);
plot(L4_axial, NRMSE_axial1(L4_axial,9)*100, '.', 'LineWidth', 1, 'Color', cmap(7-4,:), 'MarkerSize', 15);
xlim([1 12]);
ylim([0 30]);
set(gca, 'Box', 'On');
xlabel('Rank L');
ylabel('Axial, NRMSE (%)');
title({'Higher-order encoding with', 'only static off-resonance'});
text(1.1, 30, {'(A)'}, 'Color', 'k', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Red box
hCircle = annotation(gcf, 'rectangle', [0.283+0.0012 0.4920+0.0021+0.004 0.0192 0.0287], 'Color',[1 1 1], 'FaceColor',[1 0 0], 'FaceAlpha', 0.5);
%[0.283+0.0004 0.4920-0.443+0.005 0.0192 0.0287]
%--------------------------------------------------------------------------
% Axial, B0=0, conc=1
%--------------------------------------------------------------------------
ax2 = subplot(2,3,2); grid on; grid minor;
hold on;
for idx = 1:length(idx_list)
    plot(1:Lmax_axial, NRMSE_axial2(:,idx_list(idx))*100, 'LineWidth', 1, 'Color', cmap(7-idx,:));
end
plot(L1_axial, NRMSE_axial2(L1_axial,9)*100, '.', 'LineWidth', 1, 'Color', cmap(7-4,:), 'MarkerSize', 15);
plot(L2_axial, NRMSE_axial2(L2_axial,9)*100, '.', 'LineWidth', 1, 'Color', cmap(7-4,:), 'MarkerSize', 15);
plot(L3_axial, NRMSE_axial2(L3_axial,9)*100, '.', 'LineWidth', 1, 'Color', cmap(7-4,:), 'MarkerSize', 15);
plot(L4_axial, NRMSE_axial2(L4_axial,9)*100, '.', 'LineWidth', 1, 'Color', cmap(7-4,:), 'MarkerSize', 15);
xlim([1 12]);
ylim([0 30]);
set(gca, 'Box', 'On');
xlabel('Rank L');
title({'Higher-order encoding with', 'only conc. fields'});
text(1.1, 30, {'(B)'}, 'Color', 'k', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Red box
hCircle = annotation(gcf, 'rectangle', [0.283+0.2605+0.0012 0.4920+0.0021+0.004 0.0192 0.0287], 'Color',[1 1 1], 'FaceColor',[1 0 0], 'FaceAlpha', 0.5);

%--------------------------------------------------------------------------
% Axial, B0=1, conc=1
%--------------------------------------------------------------------------
ax3 = subplot(2,3,3); grid on; grid minor;
hold on;
for idx = 1:length(idx_list)
    plot(1:Lmax_axial, NRMSE_axial3(:,idx_list(idx))*100, 'LineWidth', 1, 'Color', cmap(7-idx,:));
end
plot(L1_axial, NRMSE_axial3(L1_axial,9)*100, '.', 'LineWidth', 1, 'Color', cmap(7-4,:), 'MarkerSize', 15);
plot(L2_axial, NRMSE_axial3(L2_axial,9)*100, '.', 'LineWidth', 1, 'Color', cmap(7-4,:), 'MarkerSize', 15);
plot(L3_axial, NRMSE_axial3(L3_axial,9)*100, '.', 'LineWidth', 1, 'Color', cmap(7-4,:), 'MarkerSize', 15);
plot(L4_axial, NRMSE_axial3(L4_axial,9)*100, '.', 'LineWidth', 1, 'Color', cmap(7-4,:), 'MarkerSize', 15);
xlim([1 12]);
ylim([0 30]);
set(gca, 'Box', 'On');
xlabel('Rank L');
title({'Higher-order encoding with', 'static off-resonance and conc. fields'});
text(1.1, 30, {'(C)'}, 'Color', 'k', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Red box
hCircle = annotation(gcf, 'rectangle', [0.283+0.2605*2+0.0005+0.0004+0.0012 0.4920+0.0021+0.004 0.0192 0.0287], 'Color',[1 1 1], 'FaceColor',[1 0 0], 'FaceAlpha', 0.5);

legend_text1 = cell(N3,1);
for idx = 1:N3
    legend_text1{idx} = sprintf('z = %5.1f mm', slice_offsets_axial(idx)*1e3);
end
hLegend1 = legend(legend_text1{idx_list});
set(hLegend1, 'Position', [0.8483 0.6000 0.12449 0.2441]);

%--------------------------------------------------------------------------
% Sagittal, B0=1, conc=0
%--------------------------------------------------------------------------
ax4 = subplot(2,3,4); grid on; grid minor;
hold on;
cmap = parula(7);
idx_list = [1; 7; 2; 8; 9; 11];
for idx = 1:length(idx_list)
    plot(1:Lmax_sagittal, NRMSE_sagittal1(:,idx_list(idx))*100, 'LineWidth', 1, 'Color', cmap(7-idx,:));
end
plot(15, NRMSE_sagittal1(15,11)*100, '.', 'LineWidth', 1, 'Color', cmap(7-6,:), 'MarkerSize', 15);
plot(20, NRMSE_sagittal1(20,11)*100, '.', 'LineWidth', 1, 'Color', cmap(7-6,:), 'MarkerSize', 15);
plot(30, NRMSE_sagittal1(30,11)*100, '.', 'LineWidth', 1, 'Color', cmap(7-6,:), 'MarkerSize', 15);
plot(50, NRMSE_sagittal1(50,11)*100, '.', 'LineWidth', 1, 'Color', cmap(7-6,:), 'MarkerSize', 15);
xlim([1 50]);
ylim([0 30]);
set(gca, 'Box', 'On');
xlabel('Rank L');
ylabel('Sagittal, NRMSE (%)');
text(1.5, 30, {'(D)'}, 'Color', 'k', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Red box
hCircle = annotation(gcf, 'rectangle', [0.283+0.0004 0.4920-0.443+0.005 0.0192 0.0287], 'Color',[1 1 1], 'FaceColor',[1 0 0], 'FaceAlpha', 0.5);

%--------------------------------------------------------------------------
% Sagittal, B0=0, conc=1
%--------------------------------------------------------------------------
ax5 = subplot(2,3,5); grid on; grid minor;
hold on;
for idx = 1:length(idx_list)
    plot(1:Lmax_sagittal, NRMSE_sagittal2(:,idx_list(idx))*100, 'LineWidth', 1, 'Color', cmap(7-idx,:));
end
plot(15, NRMSE_sagittal2(15,11)*100, '.', 'LineWidth', 1, 'Color', cmap(7-6,:), 'MarkerSize', 15);
plot(20, NRMSE_sagittal2(20,11)*100, '.', 'LineWidth', 1, 'Color', cmap(7-6,:), 'MarkerSize', 15);
plot(30, NRMSE_sagittal2(30,11)*100, '.', 'LineWidth', 1, 'Color', cmap(7-6,:), 'MarkerSize', 15);
plot(50, NRMSE_sagittal2(50,11)*100, '.', 'LineWidth', 1, 'Color', cmap(7-6,:), 'MarkerSize', 15);
xlim([1 50]);
ylim([0 30]);
set(gca, 'Box', 'On');
xlabel('Rank L');
text(1.5, 30, {'(E)'}, 'Color', 'k', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Red box
hCircle = annotation(gcf, 'rectangle', [0.283+0.2605+0.0004*2 0.4920-0.443+0.005 0.0192 0.0287], 'Color',[1 1 1], 'FaceColor',[1 0 0], 'FaceAlpha', 0.5);

%--------------------------------------------------------------------------
% Sagittal, B0=1, conc=1
%--------------------------------------------------------------------------
ax6 = subplot(2,3,6); grid on; grid minor;
hold on;
cmap = parula(7);
idx_list = [1; 7; 2; 8; 9; 11];
for idx = 1:length(idx_list)
    plot(1:Lmax_sagittal, NRMSE_sagittal3(:,idx_list(idx))*100, 'LineWidth', 1, 'Color', cmap(7-idx,:));
end
plot(15, NRMSE_sagittal3(15,11)*100, '.', 'LineWidth', 1, 'Color', cmap(7-6,:), 'MarkerSize', 15);
plot(20, NRMSE_sagittal3(20,11)*100, '.', 'LineWidth', 1, 'Color', cmap(7-6,:), 'MarkerSize', 15);
plot(30, NRMSE_sagittal3(30,11)*100, '.', 'LineWidth', 1, 'Color', cmap(7-6,:), 'MarkerSize', 15);
plot(50, NRMSE_sagittal3(50,11)*100, '.', 'LineWidth', 1, 'Color', cmap(7-6,:), 'MarkerSize', 15);
xlim([1 50]);
ylim([0 30]);
set(gca, 'Box', 'On');
xlabel('Rank L');
text(1.5, 30, {'(F)'}, 'Color', 'k', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Red box
hCircle = annotation(gcf, 'rectangle', [0.283+0.2605*2+0.0004*2 0.4920-0.443+0.005 0.0192 0.0287], 'Color',[1 1 1], 'FaceColor',[1 0 0], 'FaceAlpha', 0.5);

legend_text2 = cell(N3,1);
for idx = 1:N3
    legend_text2{idx} = sprintf('x = %5.1f mm', slice_offsets_sagittal(idx)*1e3);
end
hLegend2 = legend(legend_text2{idx_list});

hLegend1.ItemTokenSize = [30,18]/5;
hLegend2.ItemTokenSize = [30,18]/5;

%set(hLegend1, 'Position', [0.8183-0.001 0.5700 0.0848 0.2900]);
%set(hLegend2, 'Position', [0.8193-0.001 0.1226 0.0848 0.2900]);
set(hLegend1, 'Position', [0.8183-0.001 0.6200 0.0848 0.17500]);
set(hLegend2, 'Position', [0.8193-0.001 0.1226+0.05 0.0848 0.1750]);

pos1 = get(ax1, 'Position');
pos2 = get(ax2, 'Position');
pos3 = get(ax3, 'Position');
pos4 = get(ax4, 'Position');
pos5 = get(ax5, 'Position');
pos6 = get(ax6, 'Position');

set(ax1, 'Position', [0.1300-0.07    0.5838-0.05    0.2134+0.02    0.3412+0.02]);
set(ax2, 'Position', [0.4108-0.09    0.5838-0.05    0.2134+0.02    0.3412+0.02]);
set(ax3, 'Position', [0.6916-0.11    0.5838-0.05    0.2134+0.02    0.3412+0.02]);
set(ax4, 'Position', [0.1300-0.07    0.1100-0.05+0.03    0.2134+0.02    0.3412+0.02]);
set(ax5, 'Position', [0.4108-0.09    0.1100-0.05+0.03    0.2134+0.02    0.3412+0.02]);
set(ax6, 'Position', [0.6916-0.11    0.1100-0.05+0.03    0.2134+0.02    0.3412+0.02]);

%--------------------------------------------------------------------------
% Inset images (axial)
%--------------------------------------------------------------------------
ax7 = axes;
imagesc(abs(diff_montage_axial1)); axis image off;
text(N2_zoom1*0, N1_zoom2*0, sprintf('L=%d', L1_axial), 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom1*1, N1_zoom2*0, sprintf('L=%d', L2_axial), 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom1*0, N1_zoom2*1, sprintf('L=%d', L3_axial), 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom1*1, N1_zoom2*1, sprintf('L=%d', L4_axial), 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
set(ax7, 'Position', [0.185-0.042 0.650-0.015 0.17 0.3412-0.11]);
caxis([0 20]);
text(N2_zoom2, 0, sprintf('Diff. image (%4.1f mm)', slice_offsets_axial(9)*1e3), 'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
colormap(gca, turbo(256));
hc1 = colorbar('southOutside');
set(hc1, 'Position', [0.1793+0.003 0.5847+0.01+0.025 0.1167-0.025 0.01365]);
set(hc1, 'FontSize', 10);

ax8 = axes;
imagesc(abs(diff_montage_axial2)); axis image off;
text(N2_zoom1*0, N1_zoom2*0, sprintf('L=%d', L1_axial), 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom1*1, N1_zoom2*0, sprintf('L=%d', L2_axial), 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom1*0, N1_zoom2*1, sprintf('L=%d', L3_axial), 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom1*1, N1_zoom2*1, sprintf('L=%d', L4_axial), 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
set(ax8, 'Position', [0.665-0.261 0.650-0.015 0.17 0.3412-0.11]);
text(N2_zoom2, 0, sprintf('Diff. image (%4.1f mm)', slice_offsets_axial(9)*1e3), 'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
colormap(gca, turbo(256));
caxis([0 20]);
hc2 = colorbar('southOutside');
set(hc2, 'Position', [0.185+0.305-0.046 0.5847+0.01+0.025 0.1167-0.025 0.01365]);
set(hc2, 'FontSize', 10);
title(hc2, '[relative scale %]', 'Position', [49.5075-12 -34.128+8 0]);

ax9 = axes;
imagesc(abs(diff_montage_axial3)); axis image off;
text(N2_zoom1*0, N1_zoom2*0, sprintf('L=%d', L1_axial), 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom1*1, N1_zoom2*0, sprintf('L=%d', L2_axial), 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom1*0, N1_zoom2*1, sprintf('L=%d', L3_axial), 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom1*1, N1_zoom2*1, sprintf('L=%d', L4_axial), 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
set(ax9 , 'Position', [0.665-0.001 0.650-0.015 0.17 0.3412-0.11]);
text(N2_zoom2, 0, sprintf('Diff. image (%4.1f mm)', slice_offsets_axial(9)*1e3), 'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
colormap(gca, turbo(256));
caxis([0 20]);
hc3 = colorbar('southOutside');
set(hc3, 'Position', [0.185+0.305+0.215 0.5847+0.01+0.025 0.1167-0.025  0.01365]);
set(hc3, 'FontSize', 10);

%--------------------------------------------------------------------------
% Inset images (sagittal)
%--------------------------------------------------------------------------
ax10 = axes;
imagesc(abs(diff_montage_sagittal1)); axis image off;
text(N2_zoom2*0, N1_zoom2*0, 'L=15', 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom2*1, N1_zoom2*0, 'L=20', 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom2*0, N1_zoom2*1, 'L=30', 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom2*1, N1_zoom2*1, 'L=50', 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
set(ax10, 'Position', [0.185-0.042 0.650-0.455-0.005 0.17 0.3412-0.11]);
caxis([0 20]);
text(N2_zoom2, 0, sprintf('Diff. image (%4.1f mm)', slice_offsets_sagittal(11)*1e3), 'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
colormap(gca, turbo(256));
hc4 = colorbar('southOutside');
set(hc4, 'Position', [0.1793+0.003 0.5847+0.01-0.435+0.02-0.005 0.1167-0.025 0.01365]);
set(hc4, 'FontSize', 10);

ax11 = axes;
imagesc(abs(diff_montage_sagittal2)); axis image off;
text(N2_zoom2*0, N1_zoom2*0, 'L=15', 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom2*1, N1_zoom2*0, 'L=20', 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom2*0, N1_zoom2*1, 'L=30', 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom2*1, N1_zoom2*1, 'L=50', 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
set(ax11, 'Position', [0.665-0.261 0.650-0.455-0.005 0.17 0.3412-0.11]);
text(N2_zoom2, 0, sprintf('Diff. image (%4.1f mm)', slice_offsets_sagittal(11)*1e3), 'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
colormap(gca, turbo(256));
caxis([0 20]);
hc5 = colorbar('southOutside');
set(hc5, 'Position', [0.185+0.305-0.046 0.5847+0.01-0.435+0.02-0.005 0.1167-0.025 0.01365]);
set(hc5, 'FontSize', 10);

ax12 = axes;
imagesc(abs(diff_montage_sagittal3)); axis image off;
text(N2_zoom2*0, N1_zoom2*0, 'L=15', 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom2*1, N1_zoom2*0, 'L=20', 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom2*0, N1_zoom2*1, 'L=30', 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
text(N2_zoom2*1, N1_zoom2*1, 'L=50', 'Color', 'w', 'FontSize', 10, 'VerticalAlignment', 'top');
set(ax12, 'Position', [0.665 0.650-0.455-0.005 0.17 0.3412-0.11]);
text(N2_zoom2, 0, sprintf('Diff. image (%4.1f mm)', slice_offsets_sagittal(11)*1e3), 'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
colormap(gca, turbo(256));
caxis([0 20]);
hc6 = colorbar('southOutside');
set(hc6, 'Position', [0.185+0.305+0.215 0.5847+0.01-0.435+0.02-0.005 0.1167-0.025 0.01365]);
set(hc6, 'FontSize', 10);

export_fig('figure2', '-r600', '-tif', '-c[180,610,90,140]'); % [top,right,bottom,left]
