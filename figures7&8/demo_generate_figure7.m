% demo_generate_figure7.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 05/22/2021, Last modified: 10/10/2021

%% Clean slate
close all; clear; clc;

%% Set the directory names
result_directory    = 'E:\lowfield_maxgirf\results_meas_MID00260_FID03643_se_spiral_1102_sag_s24';
nufft_fullpath      = fullfile(result_directory, 'meas_MID00260_FID03643_se_spiral_1102_sag_s24_nufft_gpu.mat');
king_fullpath       = fullfile(result_directory, 'meas_MID00260_FID03643_se_spiral_1102_sag_s24_king.mat');
sense_fullpath      = fullfile(result_directory, 'meas_MID00260_FID03643_se_spiral_1102_sag_s24_sense_gpu.mat');
maxgirf_cp_fullpath = fullfile(result_directory, 'meas_MID00260_FID03643_se_spiral_1102_sag_s24_cpr_gpu.mat');
maxgirf_cg_fullpath = fullfile(result_directory, 'meas_MID00260_FID03643_se_spiral_1102_sag_s24_maxgirf_gpu_supp1_iter45');
cartesian_fullpath  = fullfile('E:\lowfield_maxgirf', 'meas_MID00258_FID03641_se_15b130_ro0');
B0map_fullpath      = 'B0map_nlinv_min1.0e-06_sagittal_ro0.mat';

%% Define stuff
slice_offsets = [-62.47249; 
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
actual_slice_nrs = [1,3,5,7,9,11,2,4,6,8,10].';

%% Load a B0 map
load(B0map_fullpath);
B0map = B0map_nlinv;

%% Define a function handle
if ~isempty(strfind(result_directory, '_sag_'))
    reorient = @(x) x;
    axis_text = 'x';
elseif ~isempty(strfind(result_directory, '_cor_'))
    reorient = @(x) x;
elseif ~isempty(strfind(result_directory, '_ax_'))
    reorient = @(x) flip(rot90(x, -1), 2);
    axis_text = 'z';
end

vec = @(x) x(:);

%% Load .mat files
slice_nr = 9;
slice_offset = slice_offsets(slice_nr); % [m]

%--------------------------------------------------------------------------
% Load .mat files
%--------------------------------------------------------------------------
load(nufft_fullpath);
load(king_fullpath);
load(maxgirf_cg_fullpath);
load(fullfile(result_directory, 'meas_MID00260_FID03643_se_spiral_1102_sag_s24_maxgirf_kp'));
load(cartesian_fullpath);

%% Calculate time-averaged concomitant fields
[N1,N2,Ns] = size(im_nufft_gpu);

k = output_maxgirf(slice_nr).k;
p = output_maxgirf(slice_nr).p;
T = output_maxgirf(slice_nr).T;
Nk = size(k,1);

%T = 2.5e-6 * Nk;
fc1 = reshape(k(end,4:end,1) * p(:,4:end).', [N1 N2]) / (2 * pi * T);

%% Zeropad in image-space
im_cartesian       = zpad(reorient(im_se(:,:,slice_nr)), [640 640]);
im_nufft           = zpad(reorient(im_nufft_gpu(:,:,slice_nr)), [640 640]);
im_king            = zpad(reorient(im_king(:,:,slice_nr)), [640 640]);
im_maxgirf_lowrank = zpad(reorient(im_maxgirf_gpu(:,:,slice_nr)), [640 640]);
B0map              = zpad(reorient(B0map(:,:,slice_nr)), [640 640]);
fc1                = zpad(fc1, [640 640]);

%% Scale images
N1 = 640;
N2 = 640;
c1 = floor(N1/2) + 1;
c2 = floor(N2/2) + 1;

reference_value = abs(im_nufft(c1,c2));

im2 = im_nufft / abs(im_nufft(c1,c2)) * reference_value / 1.2;
im3 = im_king / abs(im_king(c1,c2)) * reference_value / 1.2;
im4 = im_maxgirf_lowrank / abs(im_maxgirf_lowrank(c1,c2)) * reference_value;

im2_max = max(abs(im2(:))); % NUFFT
im3_max = max(abs(im3(:))); % King
im4_max = max(abs(im4(:))); % MaxGIRF
max_all = max(cat(1, im2_max, im3_max, im4_max));
cmax = max_all * 0.55;
im1 = im_cartesian / max(abs(im_cartesian(:))) * max_all;

N1_zoom = 400;
N2_zoom = 320;
offset = 40;
cmax_cartesian = 0.5;

idx1_range = (-floor(N1_zoom/2):ceil(N1_zoom/2)-1).' + floor(N1/2) + 1 + offset;
idx2_range = (-floor(N2_zoom/2):ceil(N2_zoom/2)-1).' + floor(N2/2) + 1;

%--------------------------------------------------------------------------
% Zoomed-in (Neck)
%--------------------------------------------------------------------------
idx1_zoom_range1 = (420:560).';
idx1_zoom_range1 = (428:552).';
idx2_zoom_range1 = (310:410).';
N1_box1 = length(idx1_zoom_range1);
N2_box1 = length(idx2_zoom_range1);

im_montage1 = cat(1, cat(2, im1(idx1_zoom_range1,idx2_zoom_range1), zeros(N1_box1,1), im2(idx1_zoom_range1,idx2_zoom_range1)), ...
                    cat(2, zeros(1,N2_box1), zeros(1,1), zeros(1,N2_box1)), ...
                    cat(2, im3(idx1_zoom_range1,idx2_zoom_range1), zeros(N1_box1,1), im4(idx1_zoom_range1,idx2_zoom_range1)));

%--------------------------------------------------------------------------
% Zommed-in (brain)
%--------------------------------------------------------------------------                
idx1_zoom_range2 = (228-35-10+10:352-35-60+10).';
idx2_zoom_range2 = (310+20+30-20:410+20-20-10).';
N1_box2 = length(idx1_zoom_range2);
N2_box2 = length(idx2_zoom_range2);

im_montage2 = cat(1, cat(2, im1(idx1_zoom_range2,idx2_zoom_range2), zeros(N1_box2,1), im2(idx1_zoom_range2,idx2_zoom_range2)), ...
                    cat(2, zeros(1,N2_box2), zeros(1,1), zeros(1,N2_box2)), ...
                    cat(2, im3(idx1_zoom_range2,idx2_zoom_range2)/1.2, zeros(N1_box2,1), im4(idx1_zoom_range2,idx2_zoom_range2)));

%% Display images
%figure('Color', 'k', 'Position', [14 64 1580 744]); % 262
%figure('Color', 'k', 'Position', [6 2 934 814]);
figure('Color', 'k', 'Position', [1619 -118 934 994]);
color_order = get(gca, 'colororder');

%--------------------------------------------------------------------------
% Cartesian reference
%--------------------------------------------------------------------------
ax1 = subplot(3,3,1); hold on;
imagesc(abs(im1(idx1_range,idx2_range))); axis image ij off;
caxis([0 cmax]);
colormap(gca, gray(256));
text(N2_zoom/2, 0, {'Cartesian reference'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(2, 0, {'(A)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Zoom-in (neck)
plot([idx2_zoom_range1(1)-idx2_range(1) idx2_zoom_range1(end)- idx2_range(1)], ...
     [idx1_zoom_range1(1)-idx1_range(1) idx1_zoom_range1(1)  - idx1_range(1)], '-', 'Color', color_order(2,:), 'LineWidth', 1);

plot([idx2_zoom_range1(1)-idx2_range(1) idx2_zoom_range1(1)   - idx2_range(1)], ...
     [idx1_zoom_range1(1)-idx1_range(1) idx1_zoom_range1(end) - idx1_range(1)], '-', 'Color', color_order(2,:), 'LineWidth', 1);

plot([idx2_zoom_range1(1)   - idx2_range(1) idx2_zoom_range1(end) - idx2_range(1)], ...
     [idx1_zoom_range1(end) - idx1_range(1) idx1_zoom_range1(end)   - idx1_range(1)], '-', 'Color', color_order(2,:), 'LineWidth', 1);

plot([idx2_zoom_range1(end) - idx2_range(1) idx2_zoom_range1(end) - idx2_range(1)], ...
     [idx1_zoom_range1(1)   - idx1_range(1) idx1_zoom_range1(end) - idx1_range(1)], '-', 'Color', color_order(2,:), 'LineWidth', 1);

% Zoom-in (brain)
plot([idx2_zoom_range2(1)-idx2_range(1) idx2_zoom_range2(end)- idx2_range(1)], ...
     [idx1_zoom_range2(1)-idx1_range(1) idx1_zoom_range2(1)  - idx1_range(1)], '-', 'Color', color_order(1,:), 'LineWidth', 1);

plot([idx2_zoom_range2(1)-idx2_range(1) idx2_zoom_range2(1)   - idx2_range(1)], ...
     [idx1_zoom_range2(1)-idx1_range(1) idx1_zoom_range2(end) - idx1_range(1)], '-', 'Color', color_order(1,:), 'LineWidth', 1);

plot([idx2_zoom_range2(1)   - idx2_range(1) idx2_zoom_range2(end) - idx2_range(1)], ...
     [idx1_zoom_range2(end) - idx1_range(1) idx1_zoom_range2(end) - idx1_range(1)], '-', 'Color', color_order(1,:), 'LineWidth', 1);

plot([idx2_zoom_range2(end) - idx2_range(1) idx2_zoom_range2(end) - idx2_range(1)], ...
     [idx1_zoom_range2(1)   - idx1_range(1) idx1_zoom_range2(end) - idx1_range(1)], '-', 'Color', color_order(1,:), 'LineWidth', 1);

%--------------------------------------------------------------------------
% NUFFT
%-------------------------------------------------------------------------- 
ax2 = subplot(3,3,2);
imagesc(abs(im2(idx1_range,idx2_range))); axis image off;
caxis([0 cmax]);
colormap(gca, gray(256));
text(N2_zoom/2, 0, {'NUFFT'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(N2_zoom/2, -35, sprintf('Spiral SE, slice at %s = %4.1f mm', axis_text, slice_offset*1e3), 'Color', color_order(2,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(2, 0, {'(B)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Zoomed-in
%--------------------------------------------------------------------------
ax3 = subplot(3,3,3);
imagesc(abs(im_montage2)); axis image;
caxis([0 cmax*1.0]);
set(ax3, 'XColor', color_order(1,:), 'YColor', color_order(1,:), 'XTick', [], 'YTick', [], 'LineWidth', 2);
colormap(gca, gray(256));
text(N2_box2, 0, {'Zoomed'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(N2_box2*1, 0, {'A'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
text(N2_box2*2, 0, {'B'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
text(N2_box2*1, N1_box2*1, {'C'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
text(N2_box2*2, N1_box2*1, {'D'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
text(2, 0, {'(E)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

hold on;
annotation(gcf, 'arrow', [0.7708 0.7601]      , [0.8501 0.8430]      , 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.7708 0.7601]-0.105, [0.8501 0.8430]      , 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.7708 0.7601]      , [0.8501 0.8430]-0.123, 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.7708 0.7601]-0.105, [0.8501 0.8430]-0.123, 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);

annotation(gcf, 'arrow', [0.6124 0.6015]      , [0.6619 0.6699]      , 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.6124 0.6015]+0.107, [0.6619 0.6699]      , 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.6124 0.6015]+0.107, [0.6619 0.6699]+0.123, 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.6124 0.6015]      , [0.6619 0.6699]+0.123, 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);

%--------------------------------------------------------------------------
% King's method without B0 correction
%--------------------------------------------------------------------------
ax4 = subplot(3,3,4);
imagesc(abs(im3(idx1_range,idx2_range))); axis image off;
caxis([0 cmax]);
colormap(gca, gray(256));
text(N2_zoom/2, 0, {'King''s method'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(2, 0, {'(C)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% MaxGIRF
%-------------------------------------------------------------------------- 
ax5 = subplot(3,3,5);
imagesc(abs(im4(idx1_range,idx2_range))); axis image off;
caxis([0 cmax]);
colormap(gca, gray(256));
text(N2_zoom/2, 0, {'MaxGIRF'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(2, 0, {'(D)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% zoomed-in
%--------------------------------------------------------------------------
ax6 = subplot(3,3,6);
imagesc(abs(im_montage1)); axis image;
caxis([0 cmax*0.6]);
set(ax6, 'XColor', color_order(2,:), 'YColor', color_order(2,:), 'XTick', [], 'YTick', [], 'LineWidth', 2);
colormap(gca, gray(256));
text(N2_box1, 0, {'Zoomed'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(N2_box1*1, 0, {'A'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
text(N2_box1*2, 0, {'B'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
text(N2_box1*1, N1_box1*1, {'C'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
text(N2_box1*2, N1_box1*1, {'D'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
text(2, 0, {'(F)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
hold on;
annotation(gcf, 'arrow', [0.6102 0.6209]      , [0.3893 0.3993]      , 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.6102 0.6209]+0.107, [0.3893 0.3993]      , 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.6102 0.6209]      , [0.3893 0.3993]+0.123, 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.6102 0.6209]+0.107, [0.3893 0.3993]+0.123, 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);

%--------------------------------------------------------------------------
% static off-resonance map
%--------------------------------------------------------------------------
ax7 = subplot(3,3,7);
imagesc(B0map(idx1_range,idx2_range)); axis image off;
caxis([-100 100]);
%caxis([-100 350]);
colormap(gca, hot(256));
text(N2_zoom/2, 0, {'Static off-resonance'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
hc1 = colorbar('WestOutside');
%set(hc1, 'FontSize', 10, 'Color', 'w', 'Position', [0.1461 0.1364 0.1536 0.0151]);
%title(hc1, '[Hz]', 'Color', 'w', 'Position', [110 0.654 0]); % [125.5275 0.654 0]
set(hc1, 'FontSize', 10, 'Color', 'w', 'Position', [0.1461-0.04 0.1138+0.02 0.0151 0.2503-0.05]);
title(hc1, '[Hz]', 'Color', 'w', 'Position', [5 138.154 0]); % [125.5275 0.654 0]
text(2, 0, {'(G)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% average concomitant fields map
%--------------------------------------------------------------------------
ax8 = subplot(3,3,8); hold on;
imagesc(fc1(idx1_range,idx2_range)); axis image ij off;
contour(gca, fc1(idx1_range,idx2_range), (0:50:350).', 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'w');
caxis([-50 350]);
colormap(gca, jet(256));
text(N2_zoom/2, 0, {'Averaged conc. fields'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
% hc2 = colorbar('EastOutside');
% set(hc2, 'FontSize', 10, 'Color', 'w', 'Position', [0.5684 0.1138+0.02 0.0151 0.2503-0.05]);
% title(hc2, '[Hz]', 'Color', 'w');
text(2, 0, {'(H)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% static off-resonance map + average concomitant fields map
%--------------------------------------------------------------------------
ax9 = subplot(3,3,9); hold on;
imagesc(B0map(idx1_range,idx2_range) + fc1(idx1_range,idx2_range)); axis image ij off;
contour(gca, B0map(idx1_range,idx2_range) + fc1(idx1_range,idx2_range), cat(1, 0, (0:50:350).'), ...
    'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'w');
caxis([-30 350]);
colormap(gca, jet(256));
text(N2_zoom/2, 0, {'Sum of (G) and (H)'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
hc2 = colorbar('EastOutside');
set(hc2, 'FontSize', 10, 'Color', 'w', 'Position', [0.7858 0.1138+0.02 0.0151 0.2503-0.05]);
title(hc2, '[Hz]', 'Color', 'w');
text(2, 0, {'(I)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

pos1 = get(ax1, 'Position');
pos2 = get(ax2, 'Position');
pos3 = get(ax3, 'Position');

pos4 = get(ax4, 'Position');
pos5 = get(ax5, 'Position');
pos6 = get(ax6, 'Position');

pos7 = get(ax7, 'Position');
pos8 = get(ax8, 'Position');
pos9 = get(ax9, 'Position');

set(ax1, 'Position', [0.1300 0.5838+0.03 0.2134 0.3412]);
set(ax2, 'Position', [0.3470 0.5838+0.03 0.2134 0.3412]);
set(ax3, 'Position', [0.5640 0.5838+0.03 0.2134 0.3412]);

set(ax4, 'Position', [0.1300 0.2090+0.132 0.2134 0.3412]);
set(ax5, 'Position', [0.3470 0.2090+0.132 0.2134 0.3412]);
set(ax6, 'Position', [0.5640 0.2090+0.132 0.2134 0.3412]);

set(ax7, 'Position', [0.1300 0.0682 0.2134 0.3412]);
set(ax8, 'Position', [0.3470 0.0682 0.2134 0.3412]);
set(ax9, 'Position', [0.5640 0.0682 0.2134 0.3412]);

export_fig('figure7', '-r400', '-tif', '-c[180,630,440,260]'); % [top,right,bottom,left]
