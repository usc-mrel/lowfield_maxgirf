% demo_generate_figure5.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 05/22/2021, Last modified: 09/19/2021

%% Clean slate
close all; clear; clc;

%% Set the directory names
output_directory = 'D:\lowfield_maxgirf\figures5&6\se_spiral_1102_ax_s24_osf1_B0correction1_concomitant_correction1_Lmax50_L8';
data_directory = 'D:\lowfield\NHLBI\data\20201102_NV_brain';
B0map_fullpath = 'D:\lowfield_maxgirf\B0map_nlinv_min1.0e-06_axial.mat';

%% Define stuff
slice_offsets = [-70; -35; 0; 35; 70; 105; -52.5; -17.5; 17.5; 52.5; 87.5] * 1e-3; % [m]
actual_slice_nrs = [1,3,5,7,9,11,2,4,6,8,10].';

%% Load cartesian images
load(fullfile(data_directory, 'd20201102_NV_brain'));
im_cartesian = flip(rot90(s24_se_15b130_tra.img,1),1); % axial

%% Load a B0 map
load(B0map_fullpath);
B0map = B0map_nlinv;

%% Define a function handle
if ~isempty(strfind(output_directory, '_sag_'))
    reorient = @(x) x;
    axis_text = 'x';
elseif ~isempty(strfind(output_directory, '_cor_'))
    reorient = @(x) x;
elseif ~isempty(strfind(output_directory, '_ax_'))
    reorient = @(x) flip(rot90(x, -1), 2);
    axis_text = 'z';
end

vec = @(x) x(:);

%% Load .mat files
slice_nr = 9;
slice_offset = slice_offsets(slice_nr); % [m]
actual_slice_nr = actual_slice_nrs(slice_nr);

%--------------------------------------------------------------------------
% Define .mat filenames
%--------------------------------------------------------------------------
nufft_filename                   = sprintf('nufft_slice%d.mat', slice_nr);
king_filename                    = sprintf('king_slice%d.mat', slice_nr);
king_with_B0_correction_filename = sprintf('king_with_B0_correction_slice%d.mat', slice_nr);
maxgirf_lowrank_filename         = sprintf('maxgirf_lowrank_slice%d.mat', slice_nr);

%--------------------------------------------------------------------------
% Load .mat files
%--------------------------------------------------------------------------
load(fullfile(output_directory, nufft_filename));
load(fullfile(output_directory, king_filename));
load(fullfile(output_directory, king_with_B0_correction_filename));
load(fullfile(output_directory, maxgirf_lowrank_filename));
[N1,N2] = size(im_nufft);

im_cartesian               = reorient(im_cartesian);
im_nufft                   = reorient(im_nufft);
im_king                    = reorient(im_king);
im_king_with_B0_correction = reorient(im_king_with_B0_correction);
im_maxgirf_lowrank         = reorient(im_maxgirf_lowrank);
B0map                      = reorient(B0map);

%% Scale images
im1 = im_cartesian(:,:,actual_slice_nr) / max(abs(vec(im_cartesian(:,:,actual_slice_nr))));
im2 = im_maxgirf_lowrank / max(abs(im_maxgirf_lowrank(:)));
im3 = im_king / max(abs(im_king(:)));
im4 = im_king_with_B0_correction / max(abs(im_king_with_B0_correction(:)));
cmax = 0.7;

N1_zoom = 280;
N2_zoom = 240;
offset = 14;
cmax_cartesian = 0.7;

idx1_range = (-floor(N1_zoom/2):ceil(N1_zoom/2)-1).' + floor(N1/2) + 1 + offset;
idx2_range = (-floor(N2_zoom/2):ceil(N2_zoom/2)-1).' + floor(N2/2) + 1;

idx1_zoom_range = (75:155).';
idx2_zoom_range = (50:110).';
N1_box = length(idx1_zoom_range);
N2_box = length(idx2_zoom_range);

im_montage = cat(1, cat(2, im1(idx1_zoom_range,idx2_zoom_range), zeros(N1_box,1), im2(idx1_zoom_range,idx2_zoom_range)), ...
                    cat(2, zeros(1,N2_box), zeros(1,1), zeros(1,N2_box)), ...
                    cat(2, im3(idx1_zoom_range,idx2_zoom_range), zeros(N1_box,1), im4(idx1_zoom_range,idx2_zoom_range)));

%% Display images
figure('Color', 'k', 'Position', [14 64 1580 744]);
color_order = get(gca, 'colororder');

%--------------------------------------------------------------------------
% Cartesian reference
%--------------------------------------------------------------------------
ax1 = subplot(2,3,1); hold on;
imagesc(abs(im1(idx1_range,idx2_range))); axis image ij off;
caxis([0 cmax_cartesian]);
colormap(gca, gray(256));
text(N2_zoom/2, 0, {'Cartesian reference'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(2, 0, {'(A)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

plot([idx2_zoom_range(1)-idx2_range(1) idx2_zoom_range(end)- idx2_range(1)], ...
     [idx1_zoom_range(1)-idx1_range(1) idx1_zoom_range(1)  - idx1_range(1)], ':', 'Color', color_order(2,:), 'LineWidth', 2);

plot([idx2_zoom_range(1)-idx2_range(1) idx2_zoom_range(1)   - idx2_range(1)], ...
     [idx1_zoom_range(1)-idx1_range(1) idx1_zoom_range(end) - idx1_range(1)], ':', 'Color', color_order(2,:), 'LineWidth', 2);

plot([idx2_zoom_range(1)   - idx2_range(1) idx2_zoom_range(end) - idx2_range(1)], ...
     [idx1_zoom_range(end) - idx1_range(1) idx1_zoom_range(end)   - idx1_range(1)], ':', 'Color', color_order(2,:), 'LineWidth', 2);

plot([idx2_zoom_range(end) - idx2_range(1) idx2_zoom_range(end) - idx2_range(1)], ...
     [idx1_zoom_range(1)   - idx1_range(1) idx1_zoom_range(end) - idx1_range(1)], ':', 'Color', color_order(2,:), 'LineWidth', 2);

%--------------------------------------------------------------------------
% MaxGIRF
%--------------------------------------------------------------------------
ax2 = subplot(2,3,2);
imagesc(abs(im2(idx1_range,idx2_range))); axis image ij off;
caxis([0 cmax]);
colormap(gca, gray(256));
text(N2_zoom/2, 0, {'MaxGIRF'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(2, 0, {'(B)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(N2_zoom/2, -25, sprintf('Spiral SE, slice at %s = %4.1f mm', axis_text, slice_offset*1e3), 'Color', color_order(2,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% Zoomed
%--------------------------------------------------------------------------
ax3 = subplot(2,3,3);
imagesc(abs(im_montage)); axis image;
caxis([0 cmax]);
set(ax3, 'XColor', color_order(2,:), 'YColor', color_order(2,:), 'XTick', [], 'YTick', [], 'LineWidth', 2);
colormap(gca, gray(256));
text(N2_box, 0, {'Zoomed'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(N2_box*1, N1_box*1, {'A'}, 'Color', color_order(1,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(N2_box*2, N1_box*1, {'B'}, 'Color', color_order(1,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(N2_box*1, N1_box*2, {'C'}, 'Color', color_order(1,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(N2_box*2, N1_box*2, {'D'}, 'Color', color_order(1,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(2, 0, {'(E)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
hold on;
annotation(gcf, 'arrow', [0.4555 0.4607]     , [0.8182 0.8037]      , 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.4555 0.4607]     , [0.8182 0.8037]-0.171, 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.4555 0.4607]+0.06, [0.8182 0.8037]-0.171, 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.4555 0.4607]+0.06, [0.8182 0.8037]      , 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);

annotation(gcf, 'arrow', [0.5449 0.5392]     , [0.7795 0.7942]      , 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.5449 0.5392]     , [0.7795 0.7942]-0.171, 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.5449 0.5392]-0.06, [0.7795 0.7942]-0.171, 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);
annotation(gcf, 'arrow', [0.5449 0.5392]-0.06, [0.7795 0.7942]      , 'Color', 'r', 'HeadSize', 6, 'HeadStyle', 'plain', 'LineWidth', 2);

%--------------------------------------------------------------------------
% King's method without B0 correction
%--------------------------------------------------------------------------
ax4 = subplot(2,3,4);
imagesc(abs(im3(idx1_range,idx2_range))); axis image off;
caxis([0 cmax]);
colormap(gca, gray(256));
text(N2_zoom/2, 0, {'King''s method w/o', 'static off-res. correction'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(2, 0, {'(C)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% King's method with B0 correction
%--------------------------------------------------------------------------
ax5 = subplot(2,3,5);
imagesc(abs(im4(idx1_range,idx2_range))); axis image off;
caxis([0 cmax]);
colormap(gca, gray(256));
text(N2_zoom/2, 0, {'King''s method w/', 'static off-res. correction'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(2, 0, {'(D)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% static off-resonance map
%--------------------------------------------------------------------------
 ax6 = subplot(2,3,6);
imagesc(B0map(idx1_range,idx2_range,actual_slice_nr)); axis image off;
caxis([-100 100]);
colormap(gca, hot(256));
text(N2_zoom/2, 0, {'Static off-resonance'}, 'Color', color_order(3,:), 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
caxis([-100 100]);
colormap(gca, hot(256));
hc = colorbar('EastOutside');
set(hc, 'FontSize', 10, 'Color', 'w', 'Position', [0.5887-0.001 0.1586+0.015 0.008 0.3014]);
title(hc, '[Hz]', 'Color', 'w');
text(2, 0, {'(F)'}, 'Color', 'w', 'FontSize', 14, 'Rotation', 0, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

pos1 = get(ax1, 'Position');
pos2 = get(ax2, 'Position');
pos3 = get(ax3, 'Position');
pos4 = get(ax4, 'Position');
pos5 = get(ax5, 'Position');
pos6 = get(ax6, 'Position');

set(ax1, 'Position', [0.1300 0.5838-0.02 0.2134 0.3412]);
set(ax2, 'Position', [0.2690 0.5838-0.02 0.2134 0.3412]);
set(ax3, 'Position', [0.4077 0.5838-0.02 0.2134 0.3412]);
set(ax4, 'Position', [0.1300 0.2090-0.05 0.2134 0.3412]);
set(ax5, 'Position', [0.2690 0.2090-0.05 0.2134 0.3412]);
set(ax6, 'Position', [0.4077 0.2090-0.05 0.2134 0.3412]);

export_fig('figure5', '-r400', '-tif', '-c[90,2520,470,1090]'); % [top,right,bottom,left]
