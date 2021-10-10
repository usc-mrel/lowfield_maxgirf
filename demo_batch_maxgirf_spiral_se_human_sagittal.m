% demo_batch_maxgirf_spiral_se_human_sagittal.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 10/01/2021, Last modified: 10/01/2021

%% Clean slate
close all; clear; clc;

%% 20201102_NV_brain, sagittal
computer_type = computer;
if strcmp(computer_type, 'PCWIN64')
    B0map_fullpath = 'D:\lowfield_maxgirf\B0map_nlinv_min1.0e-06_sagittal.mat';
    data_directory = 'D:\lowfield\NHLBI\data\20201102_NV_brain';
elseif strcmp(computer_type, 'GLNXA64')
    B0map_fullpath = '/server/home/nlee/lowfield_maxgirf/B0map_nlinv_min1.0e-06_sagittal.mat';
    data_directory = '/scratch/nlee/NHLBI/data/20201102_NV_brain';
end

%% Set reconstruction parameters
osf     = 2;   % grid oversampling factor
maxiter = 15;  % number of CG iterations
Lmax    = 50;  % maximum rank of the SVD approximation of a higher-order encoding matrix
L       = 30;  % rank of the SVD approximation of a higher-order encoding matrix

%% Define parameters
spiral_filename = 'meas_MID00260_FID03643_se_spiral_1102_sag_s24';
output_directory_filename = 'figures7&8';
user_opts.vds_factor   = 75;
user_opts.discard_pre  = 20;
user_opts.discard_post = 20;
support_constraint = 1;
channel_range = [(1:14).';16;17];

% Very important parameters when you don't have a twix dataset
% These parameters are a function of a scan prescription and only applies
% to this particular dataset
main_orientation = 0; % SAGITTAL
PE_sign = 1;
RO_sign = -1;

%% With B0 correction and without concomitant field correction
static_B0_correction = 1;
concomitant_field_correction = 0;
demo_maxgirf_spiral_se_human;

%% Without B0 correction and with concomitant field correction
static_B0_correction = 0;
concomitant_field_correction = 1;
demo_maxgirf_spiral_se_human;

%% With B0 correction and with concomitant field correction
static_B0_correction = 1;
concomitant_field_correction = 1;
demo_maxgirf_spiral_se_human;
