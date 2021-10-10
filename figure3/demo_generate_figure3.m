% demo_generate_figure2.m
% Written by Namgyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 09/17/2020, Last modified: 09/17/2020

%% Clean slate
close all; clear all; clc;

%% Set directory names
%--------------------------------------------------------------------------
% Package directory
%--------------------------------------------------------------------------
src_directory = 'D:\lowfield_maxgirf_internal';

%% Add paths
addpath(genpath(src_directory));

%% Set directory names
data_directory_B01 = 'D:\lowfield_maxgirf\figure2\FISP_sagittal_s135_osf1_res0.94_SNRInf_0.55T_G24_S144.14_T9.2_Ni20_offset0_L50_Lmax80';
data_directory_B02 = 'D:\lowfield_maxgirf\figure2\FISP_sagittal_s135_osf1_res0.94_SNRInf_1.5T_G24_S144.14_T9.2_Ni20_offset0_L50_Lmax80';
data_directory_B03 = 'D:\lowfield_maxgirf\figure2\FISP_sagittal_s135_osf1_res0.94_SNRInf_3T_G24_S144.14_T9.2_Ni20_offset0_L50_Lmax80';
data_directory_B04 = 'D:\lowfield_maxgirf\figure2\FISP_sagittal_s135_osf1_res0.94_SNRInf_7T_G24_S144.14_T9.2_Ni20_offset0_L50_Lmax80';

data_directory_iso1 = 'D:\lowfield_maxgirf\figure2\FISP_sagittal_s135_osf1_res0.94_SNRInf_0.55T_G24_S144.14_T9.2_Ni20_offset0_L50_Lmax80';
data_directory_iso2 = 'D:\lowfield_maxgirf\figure2\FISP_sagittal_s135_osf1_res0.94_SNRInf_0.55T_G24_S144.14_T9.2_Ni20_offset50_L50_Lmax80';
data_directory_iso3 = 'D:\lowfield_maxgirf\figure2\FISP_sagittal_s135_osf1_res0.94_SNRInf_0.55T_G24_S144.14_T9.2_Ni20_offset100_L50_Lmax80';
data_directory_iso4 = 'D:\lowfield_maxgirf\figure2\FISP_sagittal_s135_osf1_res0.94_SNRInf_0.55T_G24_S144.14_T9.2_Ni20_offset150_L50_Lmax80';
data_directory_iso5 = 'D:\lowfield_maxgirf\figure2\FISP_sagittal_s135_osf1_res0.94_SNRInf_0.55T_G24_S144.14_T9.2_Ni20_offset200_L50_Lmax80';

%% Read .mat files for B0 depedence (0.55T)
%--------------------------------------------------------------------------
% Ground truth
%--------------------------------------------------------------------------
start_time = tic;
fprintf('Reading %s... ', fullfile(data_directory_B01, 'ground_truth'));
load(fullfile(data_directory_B01, 'ground_truth'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_ground_truth = flip(im_ground_truth,1);

%--------------------------------------------------------------------------
% KP
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B01, 'KP'));
load(fullfile(data_directory_B01, 'KP'));
fprintf('done! (%5.3f sec)\n', toc(start_time));

%--------------------------------------------------------------------------
% NUFFT
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B01, 'nufft'));
load(fullfile(data_directory_B01, 'nufft'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_nufft1 = flip(im_nufft,1);

%--------------------------------------------------------------------------
% King's method
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B01, 'king'));
load(fullfile(data_directory_B01, 'king'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_king1 = flip(im_king,1);

%--------------------------------------------------------------------------
% CG-based MaxGIRF
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B01, 'maxgirf_lowrank'));
load(fullfile(data_directory_B01, 'maxgirf_lowrank'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_maxgirf1 = flip(im_maxgirf,1);

%% Calculate time-averaged concomitant fields
[N1,N2] = size(im_ground_truth);
Nk = size(k,1);
T = 2.5e-6 * Nk;
fc1_B0 = flip(reshape(k(end,4:end,1) * p_recon(:,4:end).', [N1 N2]) / (2 * pi * T), 1);

%% Read .mat files for B0 depedence (1.5T)
%--------------------------------------------------------------------------
% Ground truth
%--------------------------------------------------------------------------
start_time = tic;
fprintf('Reading %s... ', fullfile(data_directory_B02, 'ground_truth'));
load(fullfile(data_directory_B02, 'ground_truth'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_ground_truth = flip(im_ground_truth,1);

%--------------------------------------------------------------------------
% KP
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B02, 'KP'));
load(fullfile(data_directory_B02, 'KP'));
fprintf('done! (%5.3f sec)\n', toc(start_time));

%--------------------------------------------------------------------------
% NUFFT
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B02, 'nufft'));
load(fullfile(data_directory_B02, 'nufft'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_nufft2 = flip(im_nufft,1);

%--------------------------------------------------------------------------
% King's method
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B02, 'king'));
load(fullfile(data_directory_B02, 'king'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_king2 = flip(im_king,1);

%--------------------------------------------------------------------------
% CG-based MaxGIRF
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B02, 'maxgirf_lowrank'));
load(fullfile(data_directory_B02, 'maxgirf_lowrank'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_maxgirf2 = flip(im_maxgirf,1);

%% Calculate time-averaged concomitant fields
[N1,N2] = size(im_ground_truth);
Nk = size(k,1);
T = 2.5e-6 * Nk;
fc2_B0 = flip(reshape(k(end,4:end,1) * p_recon(:,4:end).', [N1 N2]) / (2 * pi * T), 1);

%% Read .mat files for B0 depedence (3T)
%--------------------------------------------------------------------------
% Ground truth
%--------------------------------------------------------------------------
start_time = tic;
fprintf('Reading %s... ', fullfile(data_directory_B03, 'ground_truth'));
load(fullfile(data_directory_B03, 'ground_truth'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_ground_truth = flip(im_ground_truth,1);

%--------------------------------------------------------------------------
% KP
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B03, 'KP'));
load(fullfile(data_directory_B03, 'KP'));
fprintf('done! (%5.3f sec)\n', toc(start_time));

%--------------------------------------------------------------------------
% NUFFT
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B03, 'nufft'));
load(fullfile(data_directory_B03, 'nufft'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_nufft3 = flip(im_nufft,1);

%--------------------------------------------------------------------------
% King's method
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B03, 'king'));
load(fullfile(data_directory_B03, 'king'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_king3 = flip(im_king,1);

%--------------------------------------------------------------------------
% CG-based MaxGIRF
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B03, 'maxgirf_lowrank'));
load(fullfile(data_directory_B03, 'maxgirf_lowrank'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_maxgirf3 = flip(im_maxgirf,1);

%% Calculate time-averaged concomitant fields
[N1,N2] = size(im_ground_truth);
Nk = size(k,1);
T = 2.5e-6 * Nk;
fc3_B0 = flip(reshape(k(end,4:end,1) * p_recon(:,4:end).', [N1 N2]) / (2 * pi * T), 1);

%% Read .mat files for B0 depedence (7T)
%--------------------------------------------------------------------------
% Ground truth
%--------------------------------------------------------------------------
start_time = tic;
fprintf('Reading %s... ', fullfile(data_directory_B04, 'ground_truth'));
load(fullfile(data_directory_B04, 'ground_truth'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_ground_truth = flip(im_ground_truth,1);

%--------------------------------------------------------------------------
% KP
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B04, 'KP'));
load(fullfile(data_directory_B04, 'KP'));
fprintf('done! (%5.3f sec)\n', toc(start_time));

%--------------------------------------------------------------------------
% NUFFT
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B04, 'nufft'));
load(fullfile(data_directory_B04, 'nufft'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_nufft4 = flip(im_nufft,1);

%--------------------------------------------------------------------------
% King's method
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B04, 'king'));
load(fullfile(data_directory_B04, 'king'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_king4 = flip(im_king,1);

%--------------------------------------------------------------------------
% CG-based MaxGIRF
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_B04, 'maxgirf_lowrank'));
load(fullfile(data_directory_B04, 'maxgirf_lowrank'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_maxgirf4 = flip(im_maxgirf,1);

%% Calculate time-averaged concomitant fields
[N1,N2] = size(im_ground_truth);
Nk = size(k,1);
T = 2.5e-6 * Nk;
fc4_B0 = flip(reshape(k(end,4:end,1) * p_recon(:,4:end).', [N1 N2]) / (2 * pi * T), 1);

%% Read .mat files for off-isocenter depedence (0 mm)
%--------------------------------------------------------------------------
% Ground truth
%--------------------------------------------------------------------------
start_time = tic;
fprintf('Reading %s... ', fullfile(data_directory_iso1, 'ground_truth'));
load(fullfile(data_directory_iso1, 'ground_truth'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_ground_truth = flip(im_ground_truth,1);

%--------------------------------------------------------------------------
% KP
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso1, 'KP'));
load(fullfile(data_directory_iso1, 'KP'));
fprintf('done! (%5.3f sec)\n', toc(start_time));

%--------------------------------------------------------------------------
% NUFFT
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso1, 'nufft'));
load(fullfile(data_directory_iso1, 'nufft'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_nufft_iso1 = flip(im_nufft,1);

%--------------------------------------------------------------------------
% King's method
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso1, 'king'));
load(fullfile(data_directory_iso1, 'king'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_king_iso1 = flip(im_king,1);

%--------------------------------------------------------------------------
% CG-based MaxGIRF
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso1, 'maxgirf_lowrank'));
load(fullfile(data_directory_iso1, 'maxgirf_lowrank'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_maxgirf_iso1 = flip(im_maxgirf,1);

%% Calculate time-averaged concomitant fields
[N1,N2] = size(im_ground_truth);
Nk = size(k,1);
T = 2.5e-6 * Nk;
fc1_iso = flip(reshape(k(end,4:end,1) * p_recon(:,4:end).', [N1 N2]) / (2 * pi * T), 1);

%% Read .mat files for off-isocenter depedence (50 mm)
%--------------------------------------------------------------------------
% Ground truth
%--------------------------------------------------------------------------
start_time = tic;
fprintf('Reading %s... ', fullfile(data_directory_iso2, 'ground_truth'));
load(fullfile(data_directory_iso2, 'ground_truth'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_ground_truth = flip(im_ground_truth,1);

%--------------------------------------------------------------------------
% KP
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso2, 'KP'));
load(fullfile(data_directory_iso2, 'KP'));
fprintf('done! (%5.3f sec)\n', toc(start_time));

%--------------------------------------------------------------------------
% NUFFT
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso2, 'nufft'));
load(fullfile(data_directory_iso2, 'nufft'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_nufft_iso2 = flip(im_nufft,1);

%--------------------------------------------------------------------------
% King's method
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso2, 'king'));
load(fullfile(data_directory_iso2, 'king'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_king_iso2 = flip(im_king,1);

%--------------------------------------------------------------------------
% CG-based MaxGIRF
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso2, 'maxgirf_lowrank'));
load(fullfile(data_directory_iso2, 'maxgirf_lowrank'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_maxgirf_iso2 = flip(im_maxgirf,1);

%% Calculate time-averaged concomitant fields
[N1,N2] = size(im_ground_truth);
Nk = size(k,1);
T = 2.5e-6 * Nk;
fc2_iso = flip(reshape(k(end,4:end,1) * p_recon(:,4:end).', [N1 N2]) / (2 * pi * T), 1);

%% Read .mat files for off-isocenter depedence (100 mm)
%--------------------------------------------------------------------------
% Ground truth
%--------------------------------------------------------------------------
start_time = tic;
fprintf('Reading %s... ', fullfile(data_directory_iso3, 'ground_truth'));
load(fullfile(data_directory_iso3, 'ground_truth'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_ground_truth = flip(im_ground_truth,1);

%--------------------------------------------------------------------------
% KP
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso3, 'KP'));
load(fullfile(data_directory_iso3, 'KP'));
fprintf('done! (%5.3f sec)\n', toc(start_time));

%--------------------------------------------------------------------------
% NUFFT
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso3, 'nufft'));
load(fullfile(data_directory_iso3, 'nufft'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_nufft_iso3 = flip(im_nufft,1);

%--------------------------------------------------------------------------
% King's method
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso3, 'king'));
load(fullfile(data_directory_iso3, 'king'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_king_iso3 = flip(im_king,1);

%--------------------------------------------------------------------------
% CG-based MaxGIRF
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso3, 'maxgirf_lowrank'));
load(fullfile(data_directory_iso3, 'maxgirf_lowrank'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_maxgirf_iso3 = flip(im_maxgirf,1);

%% Calculate time-averaged concomitant fields
[N1,N2] = size(im_ground_truth);
Nk = size(k,1);
T = 2.5e-6 * Nk;
fc3_iso = flip(reshape(k(end,4:end,1) * p_recon(:,4:end).', [N1 N2]) / (2 * pi * T), 1);

%% Read .mat files for off-isocenter depedence (150 mm)
%--------------------------------------------------------------------------
% Ground truth
%--------------------------------------------------------------------------
start_time = tic;
fprintf('Reading %s... ', fullfile(data_directory_iso4, 'ground_truth'));
load(fullfile(data_directory_iso4, 'ground_truth'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_ground_truth = flip(im_ground_truth,1);

%--------------------------------------------------------------------------
% KP
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso4, 'KP'));
load(fullfile(data_directory_iso4, 'KP'));
fprintf('done! (%5.3f sec)\n', toc(start_time));

%--------------------------------------------------------------------------
% NUFFT
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso4, 'nufft'));
load(fullfile(data_directory_iso4, 'nufft'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_nufft_iso4 = flip(im_nufft,1);

%--------------------------------------------------------------------------
% King's method
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso4, 'king'));
load(fullfile(data_directory_iso4, 'king'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_king_iso4 = flip(im_king,1);

%--------------------------------------------------------------------------
% CG-based MaxGIRF
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso4, 'maxgirf_lowrank'));
load(fullfile(data_directory_iso4, 'maxgirf_lowrank'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_maxgirf_iso4 = flip(im_maxgirf,1);

%% Calculate time-averaged concomitant fields
[N1,N2] = size(im_ground_truth);
Nk = size(k,1);
T = 2.5e-6 * Nk;
fc4_iso = flip(reshape(k(end,4:end,1) * p_recon(:,4:end).', [N1 N2]) / (2 * pi * T), 1);

%% Read .mat files for off-isocenter depedence (200 mm)
%--------------------------------------------------------------------------
% Ground truth
%--------------------------------------------------------------------------
start_time = tic;
fprintf('Reading %s... ', fullfile(data_directory_iso5, 'ground_truth'));
load(fullfile(data_directory_iso5, 'ground_truth'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_ground_truth = flip(im_ground_truth,1);

%--------------------------------------------------------------------------
% KP
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso5, 'KP'));
load(fullfile(data_directory_iso5, 'KP'));
fprintf('done! (%5.3f sec)\n', toc(start_time));

%--------------------------------------------------------------------------
% NUFFT
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso5, 'nufft'));
load(fullfile(data_directory_iso5, 'nufft'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_nufft_iso5 = flip(im_nufft,1);

%--------------------------------------------------------------------------
% King's method
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso5, 'king'));
load(fullfile(data_directory_iso5, 'king'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_king_iso5 = flip(im_king,1);

%--------------------------------------------------------------------------
% CG-based MaxGIRF
%--------------------------------------------------------------------------
fprintf('Reading %s... ', fullfile(data_directory_iso5, 'maxgirf_lowrank'));
load(fullfile(data_directory_iso5, 'maxgirf_lowrank'));
fprintf('done! (%5.3f sec)\n', toc(start_time));
im_maxgirf_iso5 = flip(im_maxgirf,1);

%% Calculate time-averaged concomitant fields
[N1,N2] = size(im_ground_truth);
Nk = size(k,1);
T = 2.5e-6 * Nk;
fc5_iso = flip(reshape(k(end,4:end,1) * p_recon(:,4:end).', [N1 N2]) / (2 * pi * T), 1);

%% Scale images
c1 = floor(N1/2) + 1;
c2 = floor(N2/2) + 1;

%--------------------------------------------------------------------------
% 0.55T
%--------------------------------------------------------------------------
im_ground_truth_scale = abs(im_ground_truth(c1,c2));
im_nufft1_scale       = abs(im_nufft1(c1,c2));
im_king1_scale        = abs(im_king1(c1,c2));
im_maxgirf1_scale     = abs(im_maxgirf1(c1,c2));

im_ground_truth_scaled = im_ground_truth / im_ground_truth_scale;
im_nufft1_scaled       = im_nufft1 / im_nufft1_scale;
im_king1_scaled        = im_king1 / im_king1_scale;
im_maxgirf1_scaled     = im_maxgirf1 / im_maxgirf1_scale;

%--------------------------------------------------------------------------
% 1.5T
%--------------------------------------------------------------------------
im_ground_truth_scale = abs(im_ground_truth(c1,c2));
im_nufft2_scale       = abs(im_nufft2(c1,c2));
im_king2_scale        = abs(im_king2(c1,c2));
im_maxgirf2_scale     = abs(im_maxgirf2(c1,c2));

im_nufft2_scaled       = im_nufft2 / im_nufft2_scale;
im_king2_scaled        = im_king2 / im_king2_scale;
im_maxgirf2_scaled     = im_maxgirf2 / im_maxgirf2_scale;

%--------------------------------------------------------------------------
% 3T
%--------------------------------------------------------------------------
im_ground_truth_scale = abs(im_ground_truth(c1,c2));
im_nufft3_scale       = abs(im_nufft3(c1,c2));
im_king3_scale        = abs(im_king3(c1,c2));
im_maxgirf3_scale     = abs(im_maxgirf3(c1,c2));

im_ground_truth_scaled = im_ground_truth / im_ground_truth_scale;
im_nufft3_scaled       = im_nufft3 / im_nufft3_scale;
im_king3_scaled        = im_king3 / im_king3_scale;
im_maxgirf3_scaled     = im_maxgirf3 / im_maxgirf3_scale;

%--------------------------------------------------------------------------
% 7T
%--------------------------------------------------------------------------
im_ground_truth_scale = abs(im_ground_truth(c1,c2));
im_nufft4_scale       = abs(im_nufft4(c1,c2));
im_king4_scale        = abs(im_king4(c1,c2));
im_maxgirf4_scale     = abs(im_maxgirf4(c1,c2));

im_ground_truth_scaled = im_ground_truth / im_ground_truth_scale;
im_nufft4_scaled       = im_nufft4 / im_nufft4_scale;
im_king4_scaled        = im_king4 / im_king4_scale;
im_maxgirf4_scaled     = im_maxgirf4 / im_maxgirf4_scale;

%--------------------------------------------------------------------------
% off-isocenter dependence (0 mm)
%--------------------------------------------------------------------------
im_ground_truth_scale = abs(im_ground_truth(c1,c2));
im_nufft_iso1_scale   = abs(im_nufft_iso1(c1,c2));
im_king_iso1_scale    = abs(im_king_iso1(c1,c2));
im_maxgirf_iso1_scale = abs(im_maxgirf_iso1(c1,c2));

im_ground_truth_scaled = im_ground_truth / im_ground_truth_scale;
im_nufft_iso1_scaled   = im_nufft_iso1 / im_nufft_iso1_scale;
im_king_iso1_scaled    = im_king_iso1 / im_king_iso1_scale;
im_maxgirf_iso1_scaled = im_maxgirf_iso1 / im_maxgirf_iso1_scale;

%--------------------------------------------------------------------------
% off-isocenter dependence (50 mm)
%--------------------------------------------------------------------------
im_ground_truth_scale = abs(im_ground_truth(c1,c2));
im_nufft_iso2_scale   = abs(im_nufft_iso2(c1,c2));
im_king_iso2_scale    = abs(im_king_iso2(c1,c2));
im_maxgirf_iso2_scale = abs(im_maxgirf_iso2(c1,c2));

im_ground_truth_scaled = im_ground_truth / im_ground_truth_scale;
im_nufft_iso2_scaled   = im_nufft_iso2 / im_nufft_iso2_scale;
im_king_iso2_scaled    = im_king_iso2 / im_king_iso2_scale;
im_maxgirf_iso2_scaled = im_maxgirf_iso2 / im_maxgirf_iso2_scale;

%--------------------------------------------------------------------------
% off-isocenter dependence (100 mm)
%--------------------------------------------------------------------------
im_ground_truth_scale = abs(im_ground_truth(c1,c2));
im_nufft_iso3_scale   = abs(im_nufft_iso3(c1,c2));
im_king_iso3_scale    = abs(im_king_iso3(c1,c2));
im_maxgirf_iso3_scale = abs(im_maxgirf_iso3(c1,c2));

im_ground_truth_scaled = im_ground_truth / im_ground_truth_scale;
im_nufft_iso3_scaled   = im_nufft_iso3 / im_nufft_iso3_scale;
im_king_iso3_scaled    = im_king_iso3 / im_king_iso3_scale;
im_maxgirf_iso3_scaled = im_maxgirf_iso3 / im_maxgirf_iso3_scale;

%--------------------------------------------------------------------------
% off-isocenter dependence (150 mm)
%--------------------------------------------------------------------------
im_ground_truth_scale = abs(im_ground_truth(c1,c2));
im_nufft_iso4_scale   = abs(im_nufft_iso4(c1,c2));
im_king_iso4_scale    = abs(im_king_iso4(c1,c2));
im_maxgirf_iso4_scale = abs(im_maxgirf_iso4(c1,c2));

im_ground_truth_scaled = im_ground_truth / im_ground_truth_scale;
im_nufft_iso4_scaled   = im_nufft_iso4 / im_nufft_iso4_scale;
im_king_iso4_scaled    = im_king_iso4 / im_king_iso4_scale;
im_maxgirf_iso4_scaled = im_maxgirf_iso4 / im_maxgirf_iso4_scale;

%--------------------------------------------------------------------------
% off-isocenter dependence (200 mm)
%--------------------------------------------------------------------------
im_ground_truth_scale = abs(im_ground_truth(c1,c2));
im_nufft_iso5_scale   = abs(im_nufft_iso5(c1,c2));
im_king_iso5_scale    = abs(im_king_iso5(c1,c2));
im_maxgirf_iso5_scale = abs(im_maxgirf_iso5(c1,c2));

im_ground_truth_scaled = im_ground_truth / im_ground_truth_scale;
im_nufft_iso5_scaled   = im_nufft_iso5 / im_nufft_iso5_scale;
im_king_iso5_scaled    = im_king_iso5 / im_king_iso5_scale;
im_maxgirf_iso5_scaled = im_maxgirf_iso5 / im_maxgirf_iso5_scale;

%% Calculate the NRMSE
nrmse_nufft_B0 = zeros(4, 1, 'double');
nrmse_nufft_B0(1) = norm(im_ground_truth_scaled(:) - im_nufft1_scaled(:)) / norm(im_ground_truth_scaled(:));
nrmse_nufft_B0(2) = norm(im_ground_truth_scaled(:) - im_nufft2_scaled(:)) / norm(im_ground_truth_scaled(:));
nrmse_nufft_B0(3) = norm(im_ground_truth_scaled(:) - im_nufft3_scaled(:)) / norm(im_ground_truth_scaled(:));
nrmse_nufft_B0(4) = norm(im_ground_truth_scaled(:) - im_nufft4_scaled(:)) / norm(im_ground_truth_scaled(:));

nrmse_maxgirf_B0 = zeros(4, 1, 'double');
nrmse_maxgirf_B0(1) = norm(im_ground_truth_scaled(:) - im_maxgirf1_scaled(:)) / norm(im_ground_truth_scaled(:));
nrmse_maxgirf_B0(2) = norm(im_ground_truth_scaled(:) - im_maxgirf2_scaled(:)) / norm(im_ground_truth_scaled(:));
nrmse_maxgirf_B0(3) = norm(im_ground_truth_scaled(:) - im_maxgirf3_scaled(:)) / norm(im_ground_truth_scaled(:));
nrmse_maxgirf_B0(4) = norm(im_ground_truth_scaled(:) - im_maxgirf4_scaled(:)) / norm(im_ground_truth_scaled(:));

nrmse_nufft_iso = zeros(5, 1, 'double');
nrmse_nufft_iso(1) = norm(im_ground_truth_scaled(:) - im_nufft_iso1_scaled(:)) / norm(im_ground_truth_scaled(:));
nrmse_nufft_iso(2) = norm(im_ground_truth_scaled(:) - im_nufft_iso2_scaled(:)) / norm(im_ground_truth_scaled(:));
nrmse_nufft_iso(3) = norm(im_ground_truth_scaled(:) - im_nufft_iso3_scaled(:)) / norm(im_ground_truth_scaled(:));
nrmse_nufft_iso(4) = norm(im_ground_truth_scaled(:) - im_nufft_iso4_scaled(:)) / norm(im_ground_truth_scaled(:));
nrmse_nufft_iso(5) = norm(im_ground_truth_scaled(:) - im_nufft_iso5_scaled(:)) / norm(im_ground_truth_scaled(:));

nrmse_maxgirf_iso = zeros(5, 1, 'double');
nrmse_maxgirf_iso(1) = norm(im_ground_truth_scaled(:) - im_maxgirf_iso1_scaled(:)) / norm(im_ground_truth_scaled(:));
nrmse_maxgirf_iso(2) = norm(im_ground_truth_scaled(:) - im_maxgirf_iso2_scaled(:)) / norm(im_ground_truth_scaled(:));
nrmse_maxgirf_iso(3) = norm(im_ground_truth_scaled(:) - im_maxgirf_iso3_scaled(:)) / norm(im_ground_truth_scaled(:));
nrmse_maxgirf_iso(4) = norm(im_ground_truth_scaled(:) - im_maxgirf_iso4_scaled(:)) / norm(im_ground_truth_scaled(:));
nrmse_maxgirf_iso(5) = norm(im_ground_truth_scaled(:) - im_maxgirf_iso5_scaled(:)) / norm(im_ground_truth_scaled(:));

%% Display images
cmax = 2;
FontSize = 14;
figure('Color', 'k', 'Position', [1605 -111 801 987]); % [17 2 1064 814]
color_order = get(gca, 'colororder');

%--------------------------------------------------------------------------
% Ground truth
%--------------------------------------------------------------------------
ax1 = subplot(6,5,1);
x = reshape(p_recon(:,1), [N1 N2]);
y = reshape(p_recon(:,2), [N1 N2]);
z = reshape(p_recon(:,3), [N1 N2]);
im_ = flip(im_ground_truth_scaled,1);

surf(reshape(x*1e3, [N1 N2]), reshape(y*1e3, [N1 N2]), reshape(z*1e3, [N1 N2]), abs(im_), 'EdgeColor', 'none');
hold on;
axis equal;
set(gca, 'zdir', 'reverse', 'FontSize', 12, 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', 'TickDir', 'out', 'TickLength', [0.0100 0.0250]*4);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
colormap(gca, gray(256));
text(0, 0, -90, {'Cartesian reference', 'x = 0 mm'}, 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-2, 'Rotation', 0);
text(0, -150, -90, {'NRMSE(%)'}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);
caxis([0 cmax]);
view(-90,0);

% Zoom-in (neck)
% trial 1
idx1_zoom_range = (1:30).' + 19;
idx2_zoom_range = (1:22).' + 110;

% trial 1
idx1_zoom_range = (1:30).' + 19;
idx2_zoom_range = (1:22).' + 110 + 30;

idx1_zoom_range_ = N1 - idx1_zoom_range;
plot3([x(idx1_zoom_range_(1),idx2_zoom_range(1)) x(idx1_zoom_range_(1),idx2_zoom_range(end))]*1e3, ... 
      [y(idx1_zoom_range_(1),idx2_zoom_range(1)) y(idx1_zoom_range_(1),idx2_zoom_range(end))]*1e3, ...
      [z(idx1_zoom_range_(1),idx2_zoom_range(1)) z(idx1_zoom_range_(1),idx2_zoom_range(end))]*1e3, '-', 'Color', color_order(2,:), 'LineWidth', 1);

plot3([x(idx1_zoom_range_(1),idx2_zoom_range(1)) x(idx1_zoom_range_(end),idx2_zoom_range(1))]*1e3, ... 
      [y(idx1_zoom_range_(1),idx2_zoom_range(1)) y(idx1_zoom_range_(end),idx2_zoom_range(1))]*1e3, ...
      [z(idx1_zoom_range_(1),idx2_zoom_range(1)) z(idx1_zoom_range_(end),idx2_zoom_range(1))]*1e3, '-', 'Color', color_order(2,:), 'LineWidth', 1);

plot3([x(idx1_zoom_range_(1),idx2_zoom_range(end)) x(idx1_zoom_range_(end),idx2_zoom_range(end))]*1e3, ... 
      [y(idx1_zoom_range_(1),idx2_zoom_range(end)) y(idx1_zoom_range_(end),idx2_zoom_range(end))]*1e3, ...
      [z(idx1_zoom_range_(1),idx2_zoom_range(end)) z(idx1_zoom_range_(end),idx2_zoom_range(end))]*1e3, '-', 'Color', color_order(2,:), 'LineWidth', 1);

plot3([x(idx1_zoom_range_(end),idx2_zoom_range(1)) x(idx1_zoom_range_(end),idx2_zoom_range(end))]*1e3, ... 
      [y(idx1_zoom_range_(end),idx2_zoom_range(1)) y(idx1_zoom_range_(end),idx2_zoom_range(end))]*1e3, ...
      [z(idx1_zoom_range_(end),idx2_zoom_range(1)) z(idx1_zoom_range_(end),idx2_zoom_range(end))]*1e3, '-', 'Color', color_order(2,:), 'LineWidth', 1);

%--------------------------------------------------------------------------
% NUFFT 0.55T
%--------------------------------------------------------------------------
ax2 = subplot(6,5,2);
imagesc(abs(im_nufft1_scaled)); axis image off;
colormap(gca, gray(256));
text(N2/2, 0, sprintf('0.55T'), 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-1);
text(0, N1/2, {'NUFFT'}, 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-2, 'Rotation', 90);
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_nufft_B0(1)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% NUFFT 1.5T
%--------------------------------------------------------------------------
ax3 = subplot(6,5,3);
imagesc(abs(im_nufft2_scaled)); axis image off;
colormap(gca, gray(256));
text(N2/2, 0, sprintf('1.5T'), 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-1);
text(N2/2, -40, sprintf('(A) B0 dependence of concomitant fields'), 'Color', 'w', 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-1, 'FontWeight', 'bold');
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_nufft_B0(2)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% NUFFT 3T
%--------------------------------------------------------------------------
ax4 = subplot(6,5,4);
imagesc(abs(im_nufft3_scaled)); axis image off;
colormap(gca, gray(256));
text(N2/2, 0, sprintf('3T'), 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-1);
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_nufft_B0(3)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% NUFFT 7T
%--------------------------------------------------------------------------
ax5 = subplot(6,5,5);
imagesc(abs(im_nufft4_scaled)); axis image off;
colormap(gca, gray(256));
text(N2/2, 0, sprintf('7T'), 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-1);
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_nufft_B0(4)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% CG-based MaxGIRF 0.55T
%--------------------------------------------------------------------------
ax7 = subplot(6,5,7);
imagesc(abs(im_maxgirf1_scaled)); axis image off;
colormap(gca, gray(256));
text(0, N1/2, {'MaxGIRF'}, 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-2, 'Rotation', 90);
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_maxgirf_B0(1)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% CG-based MaxGIRF 1.5T
%--------------------------------------------------------------------------
ax8 = subplot(6,5,8);
imagesc(abs(im_maxgirf2_scaled)); axis image off;
colormap(gca, gray(256));
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_maxgirf_B0(2)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% CG-based MaxGIRF 3T
%--------------------------------------------------------------------------
ax9 = subplot(6,5,9);
imagesc(abs(im_maxgirf3_scaled)); axis image off;
colormap(gca, gray(256));
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_maxgirf_B0(3)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% CG-based MaxGIRF 7T
%--------------------------------------------------------------------------
ax10 = subplot(6,5,10);
imagesc(abs(im_maxgirf4_scaled)); axis image off;
colormap(gca, gray(256));
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_maxgirf_B0(4)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% Time-averaged concomitant fields 0.55T
%--------------------------------------------------------------------------
ax12 = subplot(6,5,12); hold on;
imagesc(abs(fc1_B0)); axis image ij off;
contour(gca, fc1_B0, cat(1, 0, (0:38:38*6).'), 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'w');
colormap(gca,jet(256));
text(0, N1/2, {'Time-averaged', 'conc. fields [Hz]'}, 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-2, 'Rotation', 90);

%--------------------------------------------------------------------------
% Time-averaged concomitant fields 1.5T
%--------------------------------------------------------------------------
ax13 = subplot(6,5,13); hold on;
imagesc(abs(fc2_B0)); axis image ij off;
contour(gca, fc2_B0, cat(1, 0, (0:14:14*6).'), 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'w');
colormap(gca,jet(256));

%--------------------------------------------------------------------------
% Time-averaged concomitant fields 3T
%--------------------------------------------------------------------------
ax14 = subplot(6,5,14); hold on;
imagesc(abs(fc3_B0)); axis image ij off;
contour(gca, fc3_B0, cat(1, 0, (0:7:7*6).'), 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'w');
colormap(gca,jet(256));

%--------------------------------------------------------------------------
% Time-averaged concomitant fields 7T
%--------------------------------------------------------------------------
ax15 = subplot(6,5,15); hold on;
imagesc(abs(fc4_B0)); axis image ij off;
contour(gca, fc4_B0, cat(1, 0, (0:3:3*6).'), 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'w');
colormap(gca,jet(256));

%--------------------------------------------------------------------------
% NUFFT (off-isocenter dependence)
%--------------------------------------------------------------------------
ax16 = subplot(6,5,16);
imagesc(abs(im_nufft_iso1_scaled)); axis image off;
colormap(gca, gray(256));
text(N2/2, 0, sprintf('x = 0 mm'), 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-1);
text(0, N1/2, {'NUFFT'}, 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-2, 'Rotation', 90);
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_nufft_iso(1)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% NUFFT (off-isocenter dependence)
%--------------------------------------------------------------------------
ax17 = subplot(6,5,17);
imagesc(abs(im_nufft_iso2_scaled)); axis image off;
colormap(gca, gray(256));
text(N2/2, 0, sprintf('x = 50 mm'), 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-1);
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_nufft_iso(2)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% NUFFT (off-isocenter dependence)
%--------------------------------------------------------------------------
ax18 = subplot(6,5,18);
imagesc(abs(im_nufft_iso3_scaled)); axis image off;
colormap(gca, gray(256));
text(N2/2, 0, sprintf('x = 100 mm'), 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-1);
text(N2/2, -40, sprintf('(B) Off-isocenter dependence of concomitant fields at 0.55T'), 'Color', 'w', 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-1, 'FontWeight', 'bold');
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_nufft_iso(3)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% NUFFT (off-isocenter dependence)
%--------------------------------------------------------------------------
ax19 = subplot(6,5,19);
imagesc(abs(im_nufft_iso4_scaled)); axis image off;
colormap(gca, gray(256));
text(N2/2, 0, sprintf('x = 150 mm'), 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-1);
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_nufft_iso(4)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% NUFFT (off-isocenter dependence)
%--------------------------------------------------------------------------
ax20 = subplot(6,5,20);
imagesc(abs(im_nufft_iso5_scaled)); axis image off;
colormap(gca, gray(256));
text(N2/2, 0, sprintf('x = 200 mm'), 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-1);
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_nufft_iso(5)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% CG-based MaxGIRF (off-isocenter dependence)
%--------------------------------------------------------------------------
ax21 = subplot(6,5,21);
imagesc(abs(im_maxgirf_iso1_scaled)); axis image off;
colormap(gca, gray(256));
text(0, N1/2, {'MaxGIRF'}, 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-2, 'Rotation', 90);
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_maxgirf_iso(1)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% CG-based MaxGIRF (off-isocenter dependence)
%--------------------------------------------------------------------------
ax22 = subplot(6,5,22);
imagesc(abs(im_maxgirf_iso2_scaled)); axis image off;
colormap(gca, gray(256));
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_maxgirf_iso(2)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% CG-based MaxGIRF (off-isocenter dependence)
%--------------------------------------------------------------------------
ax23 = subplot(6,5,23);
imagesc(abs(im_maxgirf_iso3_scaled)); axis image off;
colormap(gca, gray(256));
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_maxgirf_iso(3)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% CG-based MaxGIRF (off-isocenter dependence)
%--------------------------------------------------------------------------
ax24 = subplot(6,5,24);
imagesc(abs(im_maxgirf_iso4_scaled)); axis image off;
colormap(gca, gray(256));
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_maxgirf_iso(4)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% CG-based MaxGIRF (off-isocenter dependence)
%--------------------------------------------------------------------------
ax25 = subplot(6,5,25);
imagesc(abs(im_maxgirf_iso5_scaled)); axis image off;
colormap(gca, gray(256));
caxis([0 cmax]);
text(N2, 0, {sprintf('%5.1f', nrmse_maxgirf_iso(5)*1e2)}, 'Color', 'g', 'Interpreter', 'tex', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize-2);

%--------------------------------------------------------------------------
% Time-averaged concomitant fields (off-isocenter dependence)
%--------------------------------------------------------------------------
ax26 = subplot(6,5,26); hold on;
imagesc(abs(fc1_iso)); axis image ij off;
contour(gca, fc1_iso, cat(1, 0, (0:38:38*6).'), 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'w');
colormap(gca,jet(256));
text(0, N1/2, {'Time-averaged', 'conc. fields [Hz]'}, 'Color', color_order(3,:), 'Interpreter', 'tex', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FontSize-2, 'Rotation', 90);

%--------------------------------------------------------------------------
% Time-averaged concomitant fields (off-isocenter dependence)
%--------------------------------------------------------------------------
ax27 = subplot(6,5,27); hold on;
imagesc(abs(fc2_iso)); axis image ij off;
contour(gca, fc2_iso, cat(1, 0, (0:38:38*6).'), 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'w');
colormap(gca,jet(256));

%--------------------------------------------------------------------------
% Time-averaged concomitant fields (off-isocenter dependence)
%--------------------------------------------------------------------------
ax28 = subplot(6,5,28); hold on;
imagesc(abs(fc3_iso)); axis image ij off;
contour(gca, fc3_iso, cat(1, 0, (0:38:38*6).'), 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'w');
colormap(gca,jet(256));

%--------------------------------------------------------------------------
% Time-averaged concomitant fields (off-isocenter dependence)
%--------------------------------------------------------------------------
ax29 = subplot(6,5,29); hold on;
imagesc(abs(fc4_iso)); axis image ij off;
contour(gca, fc4_iso, cat(1, 0, (0:38:38*6).'), 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'w');
colormap(gca,jet(256));

%--------------------------------------------------------------------------
% Time-averaged concomitant fields (off-isocenter dependence)
%--------------------------------------------------------------------------
ax30 = subplot(6,5,30); hold on;
imagesc(abs(fc4_iso)); axis image ij off;
contour(gca, fc5_iso, cat(1, 0, (0:38:38*6).'), 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'w');
colormap(gca,jet(256));

pos1 = get(ax1, 'Position');
pos2 = get(ax2, 'Position');
pos3 = get(ax3, 'Position');
pos4 = get(ax4, 'Position');
pos5 = get(ax5, 'Position');

pos7 = get(ax7, 'Position');
pos8 = get(ax8, 'Position');
pos9 = get(ax9, 'Position');
pos10 = get(ax10, 'Position');

pos12 = get(ax12, 'Position');
pos13 = get(ax13, 'Position');
pos14 = get(ax14, 'Position');
pos15 = get(ax15, 'Position');

pos16 = get(ax16, 'Position');
pos17 = get(ax17, 'Position');
pos18 = get(ax18, 'Position');
pos19 = get(ax19, 'Position');
pos20 = get(ax20, 'Position');

pos21 = get(ax21, 'Position');
pos22 = get(ax22, 'Position');
pos23 = get(ax23, 'Position');
pos24 = get(ax24, 'Position');
pos25 = get(ax25, 'Position');

pos26 = get(ax26, 'Position');
pos27 = get(ax27, 'Position');
pos28 = get(ax28, 'Position');
pos29 = get(ax29, 'Position');
pos30 = get(ax30, 'Position');

set(ax1, 'Position', [0.1300-0.04       0.8224-0.017    0.1237+0.05    0.1026+0.05]);
set(ax2, 'Position', [0.2928+0.013*0    0.8224-0.017    0.1237+0.05    0.1026+0.05]);
set(ax3, 'Position', [0.4556+0.013*1    0.8224-0.017    0.1237+0.05    0.1026+0.05]);
set(ax4, 'Position', [0.6184+0.013*2    0.8224-0.017    0.1237+0.05    0.1026+0.05]);
set(ax5, 'Position', [0.7813+0.013*3    0.8224-0.017    0.1237+0.05    0.1026+0.05]);

set(ax7 , 'Position', [0.2928+0.013*0    0.6799-0.017    0.1237+0.05    0.1026+0.05]);
set(ax8 , 'Position', [0.4556+0.013*1    0.6799-0.017    0.1237+0.05    0.1026+0.05]);
set(ax9 , 'Position', [0.6184+0.013*2    0.6799-0.017    0.1237+0.05    0.1026+0.05]);
set(ax10, 'Position', [0.7813+0.013*3    0.6799-0.017    0.1237+0.05    0.1026+0.05]);

set(ax12, 'Position', [0.2928+0.013*0    0.5374-0.017    0.1237+0.05    0.1026+0.05]);
set(ax13, 'Position', [0.4556+0.013*1    0.5374-0.017    0.1237+0.05    0.1026+0.05]);
set(ax14, 'Position', [0.6184+0.013*2    0.5374-0.017    0.1237+0.05    0.1026+0.05]);
set(ax15, 'Position', [0.7813+0.013*3    0.5374-0.017    0.1237+0.05    0.1026+0.05]);

set(ax16 , 'Position', [0.1300+0.013*-1    0.3950-0.08   0.1237+0.05    0.1026+0.05]);
set(ax17 , 'Position', [0.2928+0.013*0     0.3950-0.08   0.1237+0.05    0.1026+0.05]);
set(ax18 , 'Position', [0.4556+0.013*1     0.3950-0.08   0.1237+0.05    0.1026+0.05]);
set(ax19 , 'Position', [0.6184+0.013*2     0.3950-0.08   0.1237+0.05    0.1026+0.05]);
set(ax20 , 'Position', [0.7813+0.013*3     0.3950-0.08   0.1237+0.05    0.1026+0.05]);

set(ax21 , 'Position', [0.1300+0.013*-1    0.2525-0.08   0.1237+0.05    0.1026+0.05]);
set(ax22 , 'Position', [0.2928+0.013*0     0.2525-0.08   0.1237+0.05    0.1026+0.05]);
set(ax23 , 'Position', [0.4556+0.013*1     0.2525-0.08   0.1237+0.05    0.1026+0.05]);
set(ax24 , 'Position', [0.6184+0.013*2     0.2525-0.08   0.1237+0.05    0.1026+0.05]);
set(ax25 , 'Position', [0.7813+0.013*3     0.2525-0.08   0.1237+0.05    0.1026+0.05]);

set(ax26 , 'Position', [0.1300+0.013*-1    0.1100-0.08   0.1237+0.05    0.1026+0.05]);
set(ax27 , 'Position', [0.2928+0.013*0     0.1100-0.08   0.1237+0.05    0.1026+0.05]);
set(ax28 , 'Position', [0.4556+0.013*1     0.1100-0.08   0.1237+0.05    0.1026+0.05]);
set(ax29 , 'Position', [0.6184+0.013*2     0.1100-0.08   0.1237+0.05    0.1026+0.05]);
set(ax30 , 'Position', [0.7813+0.013*3     0.1100-0.08   0.1237+0.05    0.1026+0.05]);

%--------------------------------------------------------------------------
% inset (B0 dependence)
%--------------------------------------------------------------------------
% Ground truth
axes('Position', [0.09 0.8893 0.0643 0.0610]); % [left bottom width height]
imagesc(abs(im_ground_truth_scaled(idx1_zoom_range,idx2_zoom_range))); axis image; colormap(gca,gray(256)); caxis([0 1.6]);
set(gca, 'XColor', color_order(2,:), 'YColor', color_order(2,:), 'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0], 'LineWidth', 1);

for idx = 1:4
    if idx == 1
        im_nufft_inset = im_nufft1_scaled;
        im_maxgirf_inset = im_maxgirf1_scaled;
    elseif idx == 2
        im_nufft_inset = im_nufft2_scaled;
        im_maxgirf_inset = im_maxgirf2_scaled;
    elseif idx == 3
        im_nufft_inset = im_nufft3_scaled;
        im_maxgirf_inset = im_maxgirf3_scaled;
    elseif idx == 4
        im_nufft_inset = im_nufft4_scaled;
        im_maxgirf_inset = im_maxgirf4_scaled;
    end
    % NUFFT
    axes('Position', [0.293+(idx-1)*0.1758 0.8893 0.0643 0.0610]); % [left bottom width height]
    imagesc(abs(im_nufft_inset(idx1_zoom_range,idx2_zoom_range))); axis image; colormap(gca, gray(256)); caxis([0 1.5]);
    set(gca, 'XColor', color_order(2,:), 'YColor', color_order(2,:), 'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0], 'LineWidth', 1);

    % MaxGIRF
    axes('Position', [0.293+(idx-1)*0.1758 0.8893-(0.1026+0.05-0.01) 0.0643 0.0610]); % [left bottom width height]
    imagesc(abs(im_maxgirf_inset(idx1_zoom_range,idx2_zoom_range))); axis image; colormap(gca, gray(256)); caxis([0 1.5]);
    set(gca, 'XColor', color_order(2,:), 'YColor', color_order(2,:), 'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0], 'LineWidth', 1);    
end

%--------------------------------------------------------------------------
% inset (off-isocenter dependence)
%--------------------------------------------------------------------------
for idx = 1:5
    if idx == 1
        im_nufft_inset = im_nufft_iso1_scaled;
        im_maxgirf_inset = im_maxgirf_iso1_scaled;
    elseif idx == 2
        im_nufft_inset = im_nufft_iso2_scaled;
        im_maxgirf_inset = im_maxgirf_iso2_scaled;
    elseif idx == 3
        im_nufft_inset = im_nufft_iso3_scaled;
        im_maxgirf_inset = im_maxgirf_iso3_scaled;
    elseif idx == 4
        im_nufft_inset = im_nufft_iso4_scaled;
        im_maxgirf_inset = im_maxgirf_iso4_scaled;
    elseif idx == 5
        im_nufft_inset = im_nufft_iso5_scaled;
        im_maxgirf_inset = im_maxgirf_iso5_scaled;
    end
    % NUFFT
    axes('Position', [0.293+(idx-1)*0.1758-0.1758 0.2893+0.1026+0.006 0.0643 0.0610]); % [left bottom width height]
    imagesc(abs(im_nufft_inset(idx1_zoom_range,idx2_zoom_range))); axis image; colormap(gca, gray(256)); caxis([0 1.5]);
    set(gca, 'XColor', color_order(2,:), 'YColor', color_order(2,:), 'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0], 'LineWidth', 1);

    % MaxGIRF
    axes('Position', [0.293+(idx-1)*0.1758-0.1758 0.2893-(0.1026-0.052-0.01)+0.006 0.0643 0.0610]); % [left bottom width height]
    imagesc(abs(im_maxgirf_inset(idx1_zoom_range,idx2_zoom_range))); axis image; colormap(gca, gray(256)); caxis([0 1.5]);
    set(gca, 'XColor', color_order(2,:), 'YColor', color_order(2,:), 'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0], 'LineWidth', 1);    
end

export_fig('figure2', '-r400', '-tif', '-c[0,0,120,30]'); % [top,right,bottom,left]);
