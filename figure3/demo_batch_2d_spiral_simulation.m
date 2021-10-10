% demo_batch_2d_spiral_simulation.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 09/27/2020, Last modified: 10/04/2021

%% Clean slate
close all; clear all; clc;

%% Set parameters for variable density spiral gradient design
Smax       = 1.4414e+04;                          % maximum slew rate [G/cm/sec] (e.g., 150 [mT/m/ms])
Gmax       = 2.4;                                 % maximum gradient [G/cm]
dt         = 2.5e-6;                              % sampling period [sec]
vds_factor = 75;                                  % variable density factor
FOV        = 24;                                  % spiral FOV
Fcoeff     = [FOV -FOV * (1 - vds_factor / 100)]; % FOV decreases linearly from 30 to 15cm
krmax     = 1 / (2 * (FOV / 256));                % [cycle/cm]
%krmax      = 1 / (2 * (FOV / 320 * 26.9 / 30));   % [cycle/cm]
res        = 1 / (2 * krmax) * 10;                % resolution [mm]
Ni         = 20;                                  % number of spiral interleaves

%% B0 dependence of concomitant fields
%--------------------------------------------------------------------------
% Set parameters
%--------------------------------------------------------------------------
Lmax       = 80;          % maximum rank of the SVD approximation of a higher-order encoding matrix
L          = 50;          % rank of the SVD approximation of a higher-order encoding matrix
image_ori  = 'sagittal';  % orientation of a slice
%B0         = 0.55;       % main field strength [T]
offset     = 0;           % offset from isocenter [mm]

B0_list = [0.55; 1.5; 3; 7];

%--------------------------------------------------------------------------
% Run simulation
%--------------------------------------------------------------------------
for ii = 1:length(B0_list)
    B0 = B0_list(ii);
    demo_2d_spiral_simulation;
end

%% off-isocenter dependence of concomitant fields at 0.55T
%--------------------------------------------------------------------------
% Set parameters
%--------------------------------------------------------------------------
Lmax       = 80;          % maximum rank of the SVD approximation of a higher-order encoding matrix
L          = 50;          % rank of the SVD approximation of a higher-order encoding matrix
image_ori  = 'sagittal';  % orientation of a slice
B0         = 0.55;        % main field strength [T]
%offset     = 0;           % offset from isocenter [mm]

offset_list = [0; 50; 100; 150; 200];

%--------------------------------------------------------------------------
% Run simulation
%--------------------------------------------------------------------------
for ii = 1:length(offset_list)
    offset = offset_list(ii);
    demo_2d_spiral_simulation;
end
