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
axial_data_directory1 = 'D:\lowfield_maxgirf\figure2\20201102_NV_brain_se_spiral_1102_ax_s24_osf1_B0correction1_concomitant_correction0_Lmax50_L8';
axial_data_directory2 = 'D:\lowfield_maxgirf\figure2\20201102_NV_brain_se_spiral_1102_ax_s24_osf1_B0correction0_concomitant_correction1_Lmax50_L8';
axial_data_directory3 = 'D:\lowfield_maxgirf\figure2\20201102_NV_brain_se_spiral_1102_ax_s24_osf1_B0correction1_concomitant_correction1_Lmax50_L8';

sagittal_data_directory1 = 'D:\lowfield_maxgirf\figure2\20201102_NV_brain_se_spiral_1102_sag_s24_osf2_B0_correction1_concomitant_correction0_Lmax50_L30';
sagittal_data_directory2 = 'D:\lowfield_maxgirf\figure2\20201102_NV_brain_se_spiral_1102_sag_s24_osf2_B0_correction0_concomitant_correction1_Lmax50_L30';
sagittal_data_directory3 = 'D:\lowfield_maxgirf\figure2\20201102_NV_brain_se_spiral_1102_sag_s24_osf2_B0_correction1_concomitant_correction1_Lmax50_L30';

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
Lmax_sagittal = 50;
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

%%
figure;
subplot(1,3,1); plot(NRMSE_sagittal1);
subplot(1,3,2); plot(NRMSE_sagittal2);
subplot(1,3,3); plot(NRMSE_sagittal3);


return

%%


%% Prepare inset axial images 
%--------------------------------------------------------------------------
% Load a .mat file
%--------------------------------------------------------------------------
maxgirf_cpr_filename = sprintf('maxgirf_cpr_slice%d.mat', 9);
load(fullfile(axial_data_directory, maxgirf_cpr_filename));

%--------------------------------------------------------------------------
% Calculate a mask defining the image and noise regions using 
% the iterative intermean algorithm
%--------------------------------------------------------------------------
m_full_axial = im_maxgirf_cpr(:,:,Lmax_axial);
[N1,N2] = size(m_full_axial);

%--------------------------------------------------------------------------
% Calculate the maximum magnitude
%--------------------------------------------------------------------------
mask_espirit = abs(m_full_axial) > 0;
im_max = abs(m_full_axial(mask_espirit));

%--------------------------------------------------------------------------
% Calculate a mask
%--------------------------------------------------------------------------
level = isodata(im_max, 'log');
mask = false(N1, N2);
mask(abs(m_full_axial) > level) = 1;


%%
% Inset image
idx1_range = (40:280).';
idx2_range = (30:320).';
N1_zoom1 = length(idx1_range);
N2_zoom1 = length(idx2_range);








phase_d1 = (angle(m_full_axial) - angle(im_maxgirf_cpr(:,:,2))) * 180 / pi .* mask;
phase_d2 = (angle(m_full_axial) - angle(im_maxgirf_cpr(:,:,4))) * 180 / pi .* mask;
phase_d3 = (angle(m_full_axial) - angle(im_maxgirf_cpr(:,:,6))) * 180 / pi .* mask;
phase_d4 = (angle(m_full_axial) - angle(im_maxgirf_cpr(:,:,8))) * 180 / pi .* mask;

phase_d1((phase_d1 > 350) | (phase_d1 < -350)) = 0;
phase_d2((phase_d2 > 350) | (phase_d2 < -350)) = 0;
phase_d3((phase_d3 > 350) | (phase_d3 < -350)) = 0;
phase_d4((phase_d4 > 350) | (phase_d4 < -350)) = 0;

nan_mask = zeros(N1,N2);
nan_mask(mask) = 1;
nan_mask(~mask) = NaN;

mag_d1 = (abs(m_full_axial) - abs(im_maxgirf_cpr(:,:,2))) ./ abs(m_full_axial) * 100 .* mask;
mag_d2 = (abs(m_full_axial) - abs(im_maxgirf_cpr(:,:,4))) ./ abs(m_full_axial) * 100 .* mask;
mag_d3 = (abs(m_full_axial) - abs(im_maxgirf_cpr(:,:,6))) ./ abs(m_full_axial) * 100 .* mask;
mag_d4 = (abs(m_full_axial) - abs(im_maxgirf_cpr(:,:,8))) ./ abs(m_full_axial) * 100 .* mask;

mag_d1 = mag_d1 .* nan_mask;
mag_d2 = mag_d2 .* nan_mask;
mag_d3 = mag_d3 .* nan_mask;
mag_d4 = mag_d4 .* nan_mask;

reorient = @(x) flip(rot90(x, -1), 2);

phase_montage_axial = cat(1, cat(2, reorient(phase_d1(idx1_range, idx2_range)), reorient(phase_d2(idx1_range, idx2_range))), ...
                       cat(2, reorient(phase_d3(idx1_range, idx2_range)), reorient(phase_d4(idx1_range, idx2_range))));
mag_montage_axial = cat(1, cat(2, reorient(mag_d1(idx1_range, idx2_range)), reorient(mag_d2(idx1_range, idx2_range))), ...
                     cat(2, reorient(mag_d3(idx1_range, idx2_range)), reorient(mag_d4(idx1_range, idx2_range))));


%%
[~,sorted_index] = sort(slice_offsets_sagittal);

% Load mat files
maxgirf_cpr_filename = sprintf('maxgirf_cpr_slice%d.mat', 11);
load(fullfile(sagittal_data_directory, maxgirf_cpr_filename));

% Calculate a mask defining the image and noise regions using the iterative intermean algorithm
m_full_sagittal = im_maxgirf_cpr(:,:,Lmax_sagittal);
[N1,N2,~] = size(m_full_sagittal);

%----------------------------------------------------------------------
% Calculate the maximum magnitude of images from all echo points
%----------------------------------------------------------------------
mask_espirit = abs(m_full_sagittal) > 0;
im_max = abs(m_full_sagittal(mask_espirit));
figure, imagesc(mask_espirit); axis image; drawnow;

%----------------------------------------------------------------------
% Calculate a mask
%----------------------------------------------------------------------
level = isodata(im_max, 'log');
mask = false(N1, N2);
mask(abs(m_full_sagittal) > level) = 1;
figure, imagesc(mask); axis image; drawnow;

%%
idx1_range = (180:450).';
idx2_range = (200:450).';
N1_zoom2 = length(idx1_range);
N2_zoom2 = length(idx2_range);

phase_d1 = (angle(m_full_sagittal) - angle(im_maxgirf_cpr(:,:,15))) * 180 / pi .* mask;
phase_d2 = (angle(m_full_sagittal) - angle(im_maxgirf_cpr(:,:,20))) * 180 / pi .* mask;
phase_d3 = (angle(m_full_sagittal) - angle(im_maxgirf_cpr(:,:,30))) * 180 / pi .* mask;
phase_d4 = (angle(m_full_sagittal) - angle(im_maxgirf_cpr(:,:,50))) * 180 / pi .* mask;

phase_d1((phase_d1 > 350) | (phase_d1 < -350)) = 0;
phase_d2((phase_d2 > 350) | (phase_d2 < -350)) = 0;
phase_d3((phase_d3 > 350) | (phase_d3 < -350)) = 0;
phase_d4((phase_d4 > 350) | (phase_d4 < -350)) = 0;

nan_mask = zeros(N1,N2);
nan_mask(mask) = 1;
nan_mask(~mask) = NaN;

mag_d1 = (abs(m_full_sagittal) - abs(im_maxgirf_cpr(:,:,15))) ./ abs(m_full_sagittal) * 100 .* mask;
mag_d2 = (abs(m_full_sagittal) - abs(im_maxgirf_cpr(:,:,20))) ./ abs(m_full_sagittal) * 100 .* mask;
mag_d3 = (abs(m_full_sagittal) - abs(im_maxgirf_cpr(:,:,30))) ./ abs(m_full_sagittal) * 100 .* mask;
mag_d4 = (abs(m_full_sagittal) - abs(im_maxgirf_cpr(:,:,50))) ./ abs(m_full_sagittal) * 100 .* mask;

mag_d1 = mag_d1 .* nan_mask;
mag_d2 = mag_d2 .* nan_mask;
mag_d3 = mag_d3 .* nan_mask;
mag_d4 = mag_d4 .* nan_mask;

reorient = @(x) x;

phase_montage_sagittal = cat(1, cat(2, reorient(phase_d1(idx1_range, idx2_range)), reorient(phase_d2(idx1_range, idx2_range))), ...
                                cat(2, reorient(phase_d3(idx1_range, idx2_range)), reorient(phase_d4(idx1_range, idx2_range))));
mag_montage_sagittal = cat(1, cat(2, reorient(mag_d1(idx1_range, idx2_range)), reorient(mag_d2(idx1_range, idx2_range))), ...
                              cat(2, reorient(mag_d3(idx1_range, idx2_range)), reorient(mag_d4(idx1_range, idx2_range))));
                          
%%
legend_text1 = cell(N3,1);
for idx = 1:N3
    legend_text1{idx} = sprintf('z = %5.1f mm', slice_offsets_axial(idx)*1e3);
end

figure('Color', 'w', 'Position', [14 15 980 758]);
color_order = get(gca, 'colororder');
ax1 = subplot(2,2,1); grid on; grid minor;
hold on;
cmap = parula(12);
for idx = 1:N3
    plot(1:Lmax_axial, NRMSEm_axial(:,sorted_index(idx))*100, 'LineWidth', 1, 'Color', cmap(12-idx,:));
end
plot(2, NRMSEm_axial(2,9)*100, '.', 'LineWidth', 1, 'Color', cmap(12-6,:), 'MarkerSize', 15);
plot(4, NRMSEm_axial(4,9)*100, '.', 'LineWidth', 1, 'Color', cmap(12-6,:), 'MarkerSize', 15);
plot(6, NRMSEm_axial(6,9)*100, '.', 'LineWidth', 1, 'Color', cmap(12-6,:), 'MarkerSize', 15);
plot(8, NRMSEm_axial(8,9)*100, '.', 'LineWidth', 1, 'Color', cmap(12-6,:), 'MarkerSize', 15);
xlim([1 12]);
ylim([0 30]);
set(gca, 'Box', 'On');
xlabel('Rank L');
ylabel('NRMSE (%)');
title('Magnitude NRMSE (Axial)');

ax2 = subplot(2,2,2); grid on; grid minor;
hold on;
for idx = 1:N3
    plot(1:Lmax_axial, NRMSEp_axial(:,sorted_index(idx))*100, 'LineWidth', 1, 'Color', cmap(12-idx,:));
end
plot(2, NRMSEp_axial(2,9)*100, '.', 'LineWidth', 1, 'Color', cmap(12-6,:), 'MarkerSize', 15);
plot(4, NRMSEp_axial(4,9)*100, '.', 'LineWidth', 1, 'Color', cmap(12-6,:), 'MarkerSize', 15);
plot(6, NRMSEp_axial(6,9)*100, '.', 'LineWidth', 1, 'Color', cmap(12-6,:), 'MarkerSize', 15);
plot(8, NRMSEp_axial(8,9)*100, '.', 'LineWidth', 1, 'Color', cmap(12-6,:), 'MarkerSize', 15);
hLegend1 = legend(legend_text1{sorted_index});
xlim([1 12]);
ylim([0 30]);
set(gca, 'Box', 'On');
xlabel('Rank L');
ylabel('NRMSE (%)');
title('Phase NRMSE (Axial)');

legend_text2 = cell(N3,1);
for idx = 1:N3
    legend_text2{idx} = sprintf('x = %5.1f mm', slice_offsets_sagittal(idx)*1e3);
end

ax3 = subplot(2,2,3); grid on; grid minor;
hold on;
cmap = parula(12);
for idx = 1:N3
    plot(1:Lmax_sagittal, NRMSEm_sagittal(:,sorted_index(idx))*100, 'LineWidth', 1, 'Color', cmap(12-idx,:));
end
plot(15, NRMSEm_sagittal(15,11)*100, '.', 'LineWidth', 1, 'Color', cmap(12-10,:), 'MarkerSize', 15);
plot(20, NRMSEm_sagittal(20,11)*100, '.', 'LineWidth', 1, 'Color', cmap(12-10,:), 'MarkerSize', 15);
plot(30, NRMSEm_sagittal(30,11)*100, '.', 'LineWidth', 1, 'Color', cmap(12-10,:), 'MarkerSize', 15);
plot(50, NRMSEm_sagittal(50,11)*100, '.', 'LineWidth', 1, 'Color', cmap(12-10,:), 'MarkerSize', 15);

xlim([1 50]);
ylim([0 30]);
set(gca, 'Box', 'On');
xlabel('Rank L');
ylabel('NRMSE (%)');
title('Magnitude NRMSE (Sagittal)');

ax4 = subplot(2,2,4); grid on; grid minor;
hold on;
for idx = 1:N3
    plot(1:Lmax_sagittal, NRMSEp_sagittal(:,sorted_index(idx))*100, 'LineWidth', 1, 'Color', cmap(12-idx,:));
end
plot(15, NRMSEp_sagittal(15,11)*100, '.', 'LineWidth', 1, 'Color', cmap(12-10,:), 'MarkerSize', 15);
plot(20, NRMSEp_sagittal(20,11)*100, '.', 'LineWidth', 1, 'Color', cmap(12-10,:), 'MarkerSize', 15);
plot(30, NRMSEp_sagittal(30,11)*100, '.', 'LineWidth', 1, 'Color', cmap(12-10,:), 'MarkerSize', 15);
plot(50, NRMSEp_sagittal(50,11)*100, '.', 'LineWidth', 1, 'Color', cmap(12-10,:), 'MarkerSize', 15);
hLegend2 = legend(legend_text2{sorted_index});
xlim([1 50]);
ylim([0 30]);
set(gca, 'Box', 'On');
xlabel('Rank L');
ylabel('NRMSE (%)');
title('Phase NRMSE (Sagittal)');

set(ax1, 'Position', [0.1538 0.5453 0.3207 0.3892]);
set(ax2, 'Position', [0.5193 0.5453 0.3207 0.3892]);
set(ax3, 'Position', [0.1538 0.0620 0.3207 0.3892]);
set(ax4, 'Position', [0.5193 0.0620 0.3207 0.3892]);

set(hLegend1, 'Position', [0.8483 0.6000 0.12449 0.2441]);
set(hLegend2, 'Position', [0.8483 0.1627 0.12449 0.2441]);

%--------------------------------------------------------------------------
% Inset images (sagittal)
%--------------------------------------------------------------------------
ax7 = axes;
imagesc(mag_montage_sagittal); axis image off;
text(N2_zoom2*0, N1_zoom2*0, 'L=15', 'Color', 'w', 'FontSize', 12, 'VerticalAlignment', 'top');
text(N2_zoom2*1, N1_zoom2*0, 'L=20', 'Color', 'w', 'FontSize', 12, 'VerticalAlignment', 'top');
text(N2_zoom2*0, N1_zoom2*1, 'L=30', 'Color', 'w', 'FontSize', 12, 'VerticalAlignment', 'top');
text(N2_zoom2*1, N1_zoom2*1, 'L=50', 'Color', 'w', 'FontSize', 12, 'VerticalAlignment', 'top');
set(ax7, 'Position', [0.30 0.1280 0.17 0.3412]);
caxis([-50 25]);
title(sprintf('Diff. image (%4.1f mm)', slice_offsets_sagittal(11)*1e3), 'FontSize', 12);
colormap(gca, parula(256));
hc1 = colorbar('southOutside');
set(hc1, 'Position',  [0.3194 0.1650 0.1378 0.0123]);
set(hc1, 'FontSize', 12);
title(hc1, '[relative scale %]', 'Position', [49.5075 -34.128 0]);

ax8 = axes;
imagesc(phase_montage_sagittal); axis image off;
text(N2_zoom2*0, N1_zoom2*0, 'L=15', 'Color', 'k', 'FontSize', 12, 'VerticalAlignment', 'top');
text(N2_zoom2*1, N1_zoom2*0, 'L=20', 'Color', 'k', 'FontSize', 12, 'VerticalAlignment', 'top');
text(N2_zoom2*0, N1_zoom2*1, 'L=30', 'Color', 'k', 'FontSize', 12, 'VerticalAlignment', 'top');
text(N2_zoom2*1, N1_zoom2*1, 'L=50', 'Color', 'k', 'FontSize', 12, 'VerticalAlignment', 'top');
set(ax8, 'Position', [0.665 0.1280 0.17 0.3412]);
title(sprintf('Diff. image (%4.1f mm)', slice_offsets_sagittal(11)*1e3), 'FontSize', 12);
colormap(gca, hsv(256));
caxis([-15 15]);
hc2 = colorbar('southOutside');
set(hc2, 'Position',  [0.665 + 0.3194 - 0.3 0.1650 0.1378 0.0123]);
set(hc2, 'FontSize', 12);
title(hc2, '[degree]', 'Position', [49.5075 -34.128 0]);

%--------------------------------------------------------------------------
% Inset images (axial)
%--------------------------------------------------------------------------
ax5 = axes;
imagesc(mag_montage_axial); axis image off;
text(N1_zoom1*0, N2_zoom1*0, 'L=2', 'Color', 'w', 'FontSize', 12, 'VerticalAlignment', 'top');
text(N1_zoom1*1, N2_zoom1*0, 'L=4', 'Color', 'w', 'FontSize', 12, 'VerticalAlignment', 'top');
text(N1_zoom1*0, N2_zoom1*1, 'L=6', 'Color', 'w', 'FontSize', 12, 'VerticalAlignment', 'top');
text(N1_zoom1*1, N2_zoom1*1, 'L=8', 'Color', 'w', 'FontSize', 12, 'VerticalAlignment', 'top');
set(ax5, 'Position', [0.30 0.596 0.17 0.3412]);
caxis([-50 25]);
title(sprintf('Diff. image (%4.1f mm)', slice_offsets_axial(9)*1e3), 'FontSize', 12);
colormap(gca, parula(256));
hc1 = colorbar('southOutside');
set(hc1, 'Position',  [0.3194 0.6180 0.1378 0.0123]);
set(hc1, 'FontSize', 12);
title(hc1, '[relative scale %]', 'Position', [49.5075 -34.128 0]);

ax6 = axes;
imagesc(phase_montage_axial); axis image off;
text(N1_zoom1*0, N2_zoom1*0, 'L=2', 'Color', 'k', 'FontSize', 12, 'VerticalAlignment', 'top');
text(N1_zoom1*1, N2_zoom1*0, 'L=4', 'Color', 'k', 'FontSize', 12, 'VerticalAlignment', 'top');
text(N1_zoom1*0, N2_zoom1*1, 'L=6', 'Color', 'k', 'FontSize', 12, 'VerticalAlignment', 'top');
text(N1_zoom1*1, N2_zoom1*1, 'L=8', 'Color', 'k', 'FontSize', 12, 'VerticalAlignment', 'top');
set(ax6, 'Position', [0.665 0.596 0.17 0.3412]);
title(sprintf('Diff. image (%4.1f mm)', slice_offsets_axial(9)*1e3), 'FontSize', 12);
colormap(gca, hsv(256));
caxis([-15 15]);
hc2 = colorbar('southOutside');
set(hc2, 'Position',  [0.665 + 0.3194 - 0.3 0.6180 0.1378 0.0123]);
set(hc2, 'FontSize', 12);
title(hc2, '[degree]', 'Position', [49.5075 -34.128 0]);

export_fig('figure_lowrank_v2', '-r600', '-tif');

