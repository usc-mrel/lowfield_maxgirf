% demo_non_cartesian_recon_human_axial.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/16/2022, Last modified: 02/01/2022

%% Clean slate
close all; clear all; clc;

%% Set source directories
computer_type = computer;
if strcmp(computer_type, 'PCWIN64')
    src_directory = 'E:\lowfield_maxgirf';
    ismrmrd_directory = 'D:\ismrmrd\ismrmrd';
elseif strcmp(computer_type, 'GLNXA64')
    src_directory = '/server/home/nlee/lowfield_maxgirf';
    ismrmrd_directory = '/server/home/nlee/ismrmrd';
end

%% Add source directories to search path
addpath(genpath(src_directory));
addpath(genpath(ismrmrd_directory));

%% Define data directory
computer_type = computer;
if strcmp(computer_type, 'PCWIN64')
    ismrmrd_noise_fullpath = 'D:\lowfield\NHLBI\data\20201102_NV_brain\noise\noise_meas_MID00275_FID03658_se_spiral_1102_ax_s24.h5';
    ismrmrd_data_fullpath  = 'D:\lowfield\NHLBI\data\20201102_NV_brain\h5\meas_MID00275_FID03658_se_spiral_1102_ax_s24.h5';
    siemens_dat_fullpath   = 'D:\lowfield\NHLBI\data\20201102_NV_brain\meas_MID00275_FID03658_se_spiral_1102_ax_s24.dat';
    B0map_fullpath         = 'E:\lowfield_maxgirf\B0map_nlinv_min1.0e-06_axial.mat';
elseif strcmp(computer_type, 'GLNXA64')
    ismrmrd_noise_fullpath = '/scratch/nlee/NHLBI/data/20201102_NV_brain/noise/noise_meas_MID00275_FID03658_se_spiral_1102_ax_s24.h5';
    ismrmrd_data_fullpath  = '/scratch/nlee/NHLBI/data/20201102_NV_brain/h5/meas_MID00275_FID03658_se_spiral_1102_ax_s24.h5';
    siemens_dat_fullpath   = '/scratch/nlee/NHLBI/data/20201102_NV_brain/meas_MID00275_FID03658_se_spiral_1102_ax_s24.dat';
    B0map_fullpath         = '/server/home/nlee/lowfield_maxgirf/B0map_nlinv_min1.0e-06_axial.mat';
end

%% Set reconstruction options
% "phase_sign" and "read_sign" can be determined only from Siemens raw data 
% format now until the ISMRMRD format includes these as part of its header
user_opts.phase_sign           =  1;
user_opts.read_sign            = -1;
user_opts.vds_factor           = 75;
user_opts.discard_pre          = 20;
user_opts.discard_post         = 20;
user_opts.N1                   = 320;              % reconstruction matrix size along the row direction
user_opts.N2                   = 320;              % reconstruction matrix size along the column direction
user_opts.max_iterations       = 45;               % maximum number of LSQR iterations
user_opts.tol                  = 1e-5;             % LSQR tolerance
user_opts.static_B0_correction = 1;                % static off-resonance correction: 1=yes, 0=no
user_opts.Lmax                 = 50;               % maximum rank of the SVD approximation of a higher-order encoding matrix
user_opts.L                    = 8;                % rank of the SVD approximation of a higher-order encoding matrix

%% Define an output filename
[filepath,filename,ext] = fileparts(ismrmrd_data_fullpath);
output_filename = sprintf('%s_%dx%d', filename, user_opts.N1, user_opts.N2);

%% Load a static off-resonance map [Hz]
load(B0map_fullpath);

if 0
%% Perform NUFFT reconstruction
[im_nufft,header_nufft,r_dcs_nufft] = siemens_gridding_recon(ismrmrd_noise_fullpath, ismrmrd_data_fullpath, siemens_dat_fullpath, user_opts);
save(sprintf('%s_nufft', output_filename), 'im_nufft', 'header_nufft', 'r_dcs_nufft', 'user_opts', '-v7.3');

%% Perform SENSE reconstruction
[im_sense,header_sense,r_dcs_sense,output_sense] = siemens_sense_recon(ismrmrd_noise_fullpath, ismrmrd_data_fullpath, siemens_dat_fullpath, user_opts);
save(sprintf('%s_sense', output_filename), 'im_sense', 'header_sense', 'r_dcs_sense', 'output_sense', 'user_opts', '-v7.3');

%% Perform CP-based MaxGIRF reconstruction
[im_cpr,header_cpr,r_dcs_cpr,output_cpr] = siemens_maxgirf_cp_recon(ismrmrd_noise_fullpath, ismrmrd_data_fullpath, siemens_dat_fullpath, B0map_nlinv, user_opts);
save(sprintf('%s_cpr', output_filename), 'im_cpr', 'header_cpr', 'r_dcs_cpr', 'output_cpr', 'user_opts', '-v7.3');

%% Perform CG-based MaxGIRF reconstruction
[im_maxgirf,header_maxgirf,r_dcs_maxgirf,output_maxgirf] = siemens_maxgirf_cg_recon(ismrmrd_noise_fullpath, ismrmrd_data_fullpath, siemens_dat_fullpath, B0map_nlinv, user_opts);
save(sprintf('%s_maxgirf', output_filename), 'im_maxgirf', 'header_maxgirf', 'r_dcs_maxgirf', 'output_maxgirf', 'user_opts', '-v7.3');

%% Perform CG-based MaxGIRF reconstruction (single-GPU)
[im_maxgirf_gpu,header_maxgirf,r_dcs_maxgirf,output_maxgirf] = siemens_maxgirf_cg_recon_single_gpu(ismrmrd_noise_fullpath, ismrmrd_data_fullpath, siemens_dat_fullpath, B0map_nlinv, user_opts);
save(sprintf('%s_maxgirf_single_gpu_supp%d_iter%d', output_filename, user_opts.support_constraint, user_opts.max_iterations), 'im_maxgirf_gpu', 'header_maxgirf', 'r_dcs_maxgirf', 'output_maxgirf', 'user_opts', '-v7.3');
end

%% Perform NUFFT reconstruction (single-GPU)
[im_nufft_gpu,header_nufft,r_dcs_nufft] = siemens_gridding_recon_gpu(ismrmrd_noise_fullpath, ismrmrd_data_fullpath, siemens_dat_fullpath, user_opts);
save(sprintf('%s_nufft_gpu', output_filename), 'im_nufft_gpu', 'header_nufft', 'r_dcs_nufft', 'user_opts', '-v7.3');

%% Perform King's method reconstruction
[im_king,header_king,r_dcs_king] = siemens_king_method_recon(ismrmrd_noise_fullpath, ismrmrd_data_fullpath, siemens_dat_fullpath, user_opts);
save(sprintf('%s_king', output_filename), 'im_king', 'header_king', 'r_dcs_king', 'user_opts', '-v7.3');

%% Perform SENSE reconstruction (single-GPU)
[im_sense_gpu,header_sense,r_dcs_sense,output_sense] = siemens_sense_recon_gpu(ismrmrd_noise_fullpath, ismrmrd_data_fullpath, siemens_dat_fullpath, user_opts);
save(sprintf('%s_sense_gpu', output_filename), 'im_sense_gpu', 'header_sense', 'r_dcs_sense', 'output_sense', 'user_opts', '-v7.3');

%% Perform CP-based MaxGIRF reconstruction (single-GPU)
[im_cpr_gpu,header_cpr,r_dcs_cpr,output_cpr] = siemens_maxgirf_cp_recon_gpu(ismrmrd_noise_fullpath, ismrmrd_data_fullpath, siemens_dat_fullpath, B0map_nlinv, user_opts);
save(sprintf('%s_cpr_gpu', output_filename), 'im_cpr_gpu', 'header_cpr', 'r_dcs_cpr', 'output_cpr', 'user_opts', '-v7.3');

if 0
%% Perform CG-based MaxGIRF reconstruction (single-GPU)
[im_maxgirf_gpu,header_maxgirf,r_dcs_maxgirf,output_maxgirf] = siemens_maxgirf_cg_recon_single_gpu(ismrmrd_noise_fullpath, ismrmrd_data_fullpath, siemens_dat_fullpath, B0map_nlinv, user_opts);
save(sprintf('%s_maxgirf_multi_gpu_supp%d_iter%d', output_filename, user_opts.support_constraint, user_opts.max_iterations), 'im_maxgirf_gpu', 'header_maxgirf', 'r_dcs_maxgirf', 'output_maxgirf', 'user_opts', '-v7.3');
end

%% Perform CG-based MaxGIRF reconstruction (multi-GPU)
[im_maxgirf_gpu,header_maxgirf,r_dcs_maxgirf,output_maxgirf] = siemens_maxgirf_cg_recon_multi_gpu(ismrmrd_noise_fullpath, ismrmrd_data_fullpath, siemens_dat_fullpath, B0map_nlinv, user_opts);
save(sprintf('%s_maxgirf_multi_gpu_supp%d_iter%d', output_filename, user_opts.support_constraint, user_opts.max_iterations), 'im_maxgirf_gpu', 'header_maxgirf', 'r_dcs_maxgirf', 'output_maxgirf', 'user_opts', '-v7.3');
