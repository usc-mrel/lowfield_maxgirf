% demo_MaxGIRF.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 05/23/2020, Last modified: 04/11/2021

%% Notation
%--------------------------------------------------------------------------
% We denote X, Y, Z as the horizontal, vertical, and through-magnet axes,
% respectively in the physical coordinate system, and the corresponding
% coordinates as x, y, and z. We also denote U, V, and W as the readout,
% phase encoding, and slice axes, respectively in the logical coordinate
% system, and the corresponding coordinates as u, v, and w.
%
%                     Anterior
%                        ^ (+y)
%                        |
%           Right        |            Left
%                        |
%                  Head  +--------> (+x)
%                       /
%                      /
%                     /
%               Foot v (+z)
%                      Posterior
%
%--------------------------------------------------------------------------


% BLAS call in linux
%https://www.mathworks.com/help/matlab/matlab_external/calling-lapack-and-blas-functions-from-mex-files.html#br26788-1

% Modify Function Name on UNIX Systems
% 
% Add an underscore character following the function name when calling LAPACK or BLAS functions on a UNIX® system. For example, to call dgemm, use:
% 
% dgemm_(arg1, arg2, ..., argn);
% 
% Or add these lines to your source code:
% 
% #if !defined(_WIN32)
% #define dgemm dgemm_
% #endif
    
    

%% Clean slate
close all; clear; clc;

%% Set directory names
%--------------------------------------------------------------------------
% Package directory
%--------------------------------------------------------------------------
src_directory = 'E:\lowfield_maxgirf';

%--------------------------------------------------------------------------
% ISMRMRD directory
%--------------------------------------------------------------------------
ismrmrd_directory = 'D:\ismrmrd\ismrmrd';

%--------------------------------------------------------------------------
% current directory
%--------------------------------------------------------------------------
current_directory = pwd;

%% Add paths
addpath(genpath(src_directory));
addpath(genpath(ismrmrd_directory));
warning off;

%% Define input filenames
%--------------------------------------------------------------------------
% axial at x = 0 [mm] and Gmax = 24 [mT/m]
%--------------------------------------------------------------------------
% data_directory = 'D:\lowfield\NHLBI\data\20200506_NIST_phantom_shift';
% gre_filenames{1} = 'meas_MID00182_FID69527_ax_0mm_ref_gre_TE1';
% gre_filenames{2} = 'meas_MID00183_FID69528_ax_0mm_ref_gre_TE2';
% gre_filenames{3} = 'meas_MID00184_FID69529_ax_0mm_ref_gre_TE3';
% gre_filenames{4} = 'meas_MID00185_FID69530_ax_0mm_ref_gre_TE4';
% gre_filenames{5} = 'meas_MID00186_FID69531_ax_0mm_ref_gre_TE5';
% gre_filenames{6} = 'meas_MID00187_FID69532_ax_0mm_ref_gre_TE6';
% spiral_filename = 'meas_MID00188_FID69533_ax_0mm_spiral_gmax24';
% user_opts.vds_factor   = 100;
% user_opts.discard_pre  = 0;
% user_opts.discard_post = 0;
% static_B0_correction = 1;
% support_constraint = 0;
% channel_range = [];
% 
% %--------------------------------------------------------------------------
% % axial at x = 0 [mm] and Gmax = 24 [mT/m]
% %--------------------------------------------------------------------------
% data_directory = 'D:\lowfield\NHLBI\data\20200506_NIST_phantom_shift';
% gre_filenames{1} = 'meas_MID00198_FID69543_ax_75mm_ref_gre_TE1';
% gre_filenames{2} = 'meas_MID00199_FID69544_ax_75mm_ref_gre_TE2';
% gre_filenames{3} = 'meas_MID00200_FID69545_ax_75mm_ref_gre_TE3';
% gre_filenames{4} = 'meas_MID00201_FID69546_ax_75mm_ref_gre_TE4';
% gre_filenames{5} = 'meas_MID00202_FID69547_ax_75mm_ref_gre_TE5';
% gre_filenames{6} = 'meas_MID00203_FID69548_ax_75mm_ref_gre_TE6';
% spiral_filename = 'meas_MID00204_FID69549_ax_75mm_spiral_gmax24';
% user_opts.vds_factor   = 100;
% user_opts.discard_pre  = 0;
% user_opts.discard_post = 0;
% static_B0_correction = 1;
% support_constraint = 0;
% channel_range = [];


%--------------------------------------------------------------------------
% 20201102_NV_brain, axial
%--------------------------------------------------------------------------
data_directory = 'D:\lowfield\NHLBI\data\20201102_NV_brain';
gre_filenames{1} = 'meas_MID00276_FID03659_gre_TE1';
gre_filenames{2} = 'meas_MID00277_FID03660_gre_TE2';
gre_filenames{3} = 'meas_MID00278_FID03661_gre_TE3';
gre_filenames{4} = 'meas_MID00279_FID03662_gre_TE4';
gre_filenames{5} = 'meas_MID00280_FID03663_gre_TE5';
spiral_filename = 'meas_MID00275_FID03658_se_spiral_1102_ax_s24';
user_opts.vds_factor   = 75;
user_opts.discard_pre  = 20;
user_opts.discard_post = 20;
static_B0_correction = 1;
support_constraint = 0;
channel_range = [];

%% Set reconstruction parameters
%--------------------------------------------------------------------------
% MaxGIRF reconstruction
%--------------------------------------------------------------------------
osf     = 1;    % grid oversampling factor
maxiter = 15;   % number of CG iterations
beta    = 1e-6; % Tikhonov regularization parameter
Lmax    = 15;   % maximum rank of the SVD approximation of a higher-order encoding matrix
L       = 15;   % rank of the SVD approximation of a higher-order encoding matrix

% NOT SURE WHY BETA HAS TO BE VERY LARGE TO TAKE AN EFFECT!!!
% 1e4 starts to work!
% need to scale the maximum of the coil images to 1!

%--------------------------------------------------------------------------
% TGV denoising
%--------------------------------------------------------------------------
user_opts_cartesian.zpad_factor = 2;
user_opts_cartesian.maxiter_tgv = 1500; % number of primal-dual iterations
user_opts_cartesian.lambda      = 1e2;  % regularization parameter

%% Define a fullpath to each filename
%--------------------------------------------------------------------------
% multi-echo GRE
%--------------------------------------------------------------------------
nr_TEs = length(gre_filenames);
gre_data_fullpaths  = cell(nr_TEs,1);
gre_noise_fullpaths = cell(nr_TEs,1);
gre_dat_fullpaths   = cell(nr_TEs,1);

for idx2 = 1:nr_TEs
    gre_data_fullpaths{idx2}  = fullfile(data_directory, 'h5'   , sprintf('%s.h5', gre_filenames{idx2}));
    gre_noise_fullpaths{idx2} = fullfile(data_directory, 'noise', sprintf('noise_%s.h5', gre_filenames{idx2}));
    gre_dat_fullpaths{idx2}   = fullfile(data_directory, sprintf('%s.dat', gre_filenames{idx2}));
end

%--------------------------------------------------------------------------
% spiral
%--------------------------------------------------------------------------
spiral_data_fullpath  = fullfile(data_directory, 'h5'   , sprintf('%s.h5', spiral_filename));
spiral_noise_fullpath = fullfile(data_directory, 'noise', sprintf('noise_%s.h5', spiral_filename));
spiral_dat_fullpath   = fullfile(data_directory, sprintf('%s.dat', spiral_filename));

%% Define an output directory
[filepath,data_filename] = fileparts(data_directory);
underscore_loc = strfind(spiral_filename, '_');
output_filename = sprintf('%s_osf%d_B0correction%d', spiral_filename(underscore_loc(3)+1:end), osf, static_B0_correction);
output_directory = fullfile(current_directory, data_filename, output_filename);
mkdir(output_directory);

%% Define constants
gamma = 4257.59 * (1e4 * 2 * pi); % gyromagnetic ratio for 1H [rad/sec/T]

SAGITTAL   = 0; % Patient axis perpendicular to the sagittal plane
CORONAL    = 1; % Patient axis perpendicular to the coronal plane
TRANSVERSE = 2; % Patient axis perpendicular to the transverse plane

%% Calculate a smoothed static off-resonance map [Hz]
if ~isempty(strfind(spiral_filename, 'sag'))
    user_opts_cartesian.main_orientation = SAGITTAL;
elseif ~isempty(strfind(spiral_filename, 'cor'))
    user_opts_cartesian.main_orientation = CORONAL;
elseif ~isempty(strfind(spiral_filename, 'ax'))
    user_opts_cartesian.main_orientation = TRANSVERSE;
end
user_opts_cartesian.output_directory = output_directory;
[im_echo,B0map,x_echo,y_echo,z_echo,TE_echo] = cartesian_B0map_recon(gre_noise_fullpaths, gre_data_fullpaths, gre_dat_fullpaths, user_opts_cartesian);

%% Read spiral k-space data
start_time = tic;
tic; fprintf('Reading an ismrmrd file: %s... ', spiral_data_fullpath);
if exist(spiral_data_fullpath, 'file')
    dset = ismrmrd.Dataset(spiral_data_fullpath, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , spiral_data_fullpath);
end

%% Get imaging parameters from the XML header
header = ismrmrd.xml.deserialize(dset.readxml);

%--------------------------------------------------------------------------
% measurement information
%--------------------------------------------------------------------------
patient_position = header.measurementInformation.patientPosition;

%--------------------------------------------------------------------------
% Sequence parameters
%--------------------------------------------------------------------------
TR         = header.sequenceParameters.TR * 1e-3;     % [msec] * [sec/1e3msec] => [sec]
TE         = header.sequenceParameters.TE * 1e-3;     % [msec] * [sec/1e3msec] => [sec]
flip_angle = header.sequenceParameters.flipAngle_deg; % [degrees]

%--------------------------------------------------------------------------
% Experimental conditions
%--------------------------------------------------------------------------
B0 = header.experimentalConditions.H1resonanceFrequency_Hz * (2 * pi / gamma); % [Hz] * [2pi rad/cycle] / [rad/sec/T] => [T]

%--------------------------------------------------------------------------
% Encoding
%--------------------------------------------------------------------------
field_of_view_mm(1) = header.encoding.encodedSpace.fieldOfView_mm.x; % RO
field_of_view_mm(2) = header.encoding.encodedSpace.fieldOfView_mm.y; % PE
field_of_view_mm(3) = header.encoding.encodedSpace.fieldOfView_mm.z; % SL
matrix_size(1)      = header.encoding.encodedSpace.matrixSize.x;     % number of samples in idx1 (r)
matrix_size(2)      = header.encoding.encodedSpace.matrixSize.y;     % number of samples in idx2 (c)
matrix_size(3)      = header.encoding.encodedSpace.matrixSize.z;     % number of samples in idx3 (s)

%--------------------------------------------------------------------------
% Trajectory description (HargreavesVDS2000)
%--------------------------------------------------------------------------
if strcmp(header.encoding.trajectory, 'spiral')
    Gmax              = header.encoding.trajectoryDescription.userParameterDouble(1).value;      % [G/cm]
    Smax              = header.encoding.trajectoryDescription.userParameterDouble(2).value;      % [G/cm/sec]
    FOV               = header.encoding.trajectoryDescription.userParameterDouble(3).value;      % [cm]
    krmax             = header.encoding.trajectoryDescription.userParameterDouble(4).value;      % [cycle/cm]
    interleaves       = header.encoding.trajectoryDescription.userParameterLong(1).value;        % interleaves
    fov_coefficients  = header.encoding.trajectoryDescription.userParameterLong(2).value;        % fov_coefficients
    sampling_time     = header.encoding.trajectoryDescription.userParameterLong(3).value * 1e-9; % SamplingTime_ns [sec]
    fprintf('Gmax          = %5.3f [G/cm]\n' , Gmax);
    fprintf('Smax          = %d [G/cm/sec]\n', Smax);
    fprintf('FOV           = %4.2f [cm]\n'   , FOV);
    fprintf('interleaves   = %d\n'           , interleaves);
    fprintf('sampling_time = %g [sec]\n'     , sampling_time);
end

%% Parse the ISMRMRD header
tic; fprintf('Parsing the ISMRMRD header... ');
raw_data = dset.readAcquisition(); % read all the acquisitions
fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

%--------------------------------------------------------------------------
% ISMRMRD header
%--------------------------------------------------------------------------
nr_samples         = double(max(raw_data.head.number_of_samples));
discard_pre        = double(max(raw_data.head.discard_pre));
discard_post       = double(max(raw_data.head.discard_post));
nr_channels        = double(max(raw_data.head.active_channels));
nr_phase_encodings = double(max(raw_data.head.idx.kspace_encode_step_1)) + 1;
nr_slice_encodings = double(max(raw_data.head.idx.kspace_encode_step_2)) + 1;
nr_averages        = double(max(raw_data.head.idx.average)) + 1;
nr_slices          = double(max(raw_data.head.idx.slice)) + 1;
nr_contrasts       = double(max(raw_data.head.idx.contrast)) + 1;
nr_phases          = double(max(raw_data.head.idx.phase)) + 1;
nr_repetitions     = double(max(raw_data.head.idx.repetition)) + 1;
nr_sets            = double(max(raw_data.head.idx.set)) + 1;
nr_segments        = double(max(raw_data.head.idx.segment)) + 1;

%--------------------------------------------------------------------------
% Get the dimensionality of the trajectory vector (0 means no trajectory)
%--------------------------------------------------------------------------
trajectory_dimensions = double(max(raw_data.head.trajectory_dimensions));

%--------------------------------------------------------------------------
% Get the dwell time in [sec]
%--------------------------------------------------------------------------
dt = single(max(raw_data.head.sample_time_us)) * 1e-6; % [usec] * [sec/1e-6 usec] => [sec]

%--------------------------------------------------------------------------
% Calculate the readout duration [sec]
%--------------------------------------------------------------------------
T = nr_samples * dt; % readout duration [sec]

fprintf('trajectory         = %s\n', header.encoding.trajectory);
fprintf('nr_samples         = %d\n', nr_samples);
fprintf('discard_pre        = %d\n', discard_pre);
fprintf('discard_post       = %d\n', discard_post);
fprintf('nr_channels        = %d\n', nr_channels);
fprintf('nr_phase_encodings = %d\n', nr_phase_encodings);
fprintf('nr_slice_encodings = %d\n', nr_slice_encodings);
fprintf('nr_averages        = %d\n', nr_averages);
fprintf('nr_slices          = %d\n', nr_slices);
fprintf('nr_contrasts       = %d\n', nr_contrasts);
fprintf('nr_phases          = %d\n', nr_phases);
fprintf('nr_repetitions     = %d\n', nr_repetitions);
fprintf('nr_sets            = %d\n', nr_sets);
fprintf('nr_segments        = %d\n', nr_segments);
fprintf('dt                 = %g [usec]\n', dt * 1e6);
fprintf('readout duration   = %g [msec]\n', T * 1e3);

%% Process noise only ismrmrd data
[Psi,inv_L] = process_noise_ismrmrd_data(spiral_noise_fullpath);

%% Read a Siemens dat file
fprintf('Reading a Siemens .dat file: %s\n', spiral_dat_fullpath);
twix = mapVBVD(spiral_dat_fullpath);
if length(twix) > 1
    twix = twix{end};
end

%% Define parameters used in reconstruction
N1 = osf * matrix_size(1); % number of samples in RCS (r)
N2 = osf * matrix_size(2); % number of samples in RCS (c)
N3 = matrix_size(3);       % number of samples in RCS (s)
Ni = nr_phase_encodings;   % number of spiral interleaves

%--------------------------------------------------------------------------
% Calculate the number of samples per spiral arm
%--------------------------------------------------------------------------
if ~isempty(user_opts.discard_pre)
    discard_pre = user_opts.discard_pre;
end

if ~isempty(user_opts.discard_post)
    discard_post = user_opts.discard_post;
end

Nk = nr_samples - discard_pre - discard_post; % number of samples in one spiral interleaf
Np = Nk * Ni;  % total number of fully-sampled k-space samples per image
N  = N1 * N2;  % total number of voxels in image-space
Nl = 19;       % number of spatial basis functions

%--------------------------------------------------------------------------
% Calculate the range of k-space indices for reconstruction
%--------------------------------------------------------------------------
start_index = discard_pre + 1;
end_index   = nr_samples - discard_post;
index_range = (start_index:end_index).';

%--------------------------------------------------------------------------
% Select the channels
%--------------------------------------------------------------------------
if isempty(channel_range)
    channel_range = (1:nr_channels).';
end
Nc = length(channel_range); % number of coils

%% Reconstruct images per slice
nr_recons = nr_slices * nr_contrasts * nr_phases * nr_repetitions * nr_sets * nr_segments;
% 9 for sagittal
%actual: 1 2 3 4 5 6 7 8 9 10 11
%slice : 1,3,5,7,9,11,2,4,6,8,10
%10 for axial

for idx2 = 10%:nr_recons
    %% Get information about the current slice
    [slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr] = ind2sub([nr_slices nr_contrasts nr_phases nr_repetitions nr_sets nr_segments], idx2);
    fprintf('(%2d/%2d): Reconstructing slice (slice = %2d, contrast = %2d, phase = %2d, repetition = %2d, set = %2d, segment = %2d)\n', idx2, nr_recons, slice_nr, contrast_nr, phase_nr, repetition_nr, set_nr, segment_nr);

    %% Get a list of profiles in the current slice (all averages)
    profile_list = find((raw_data.head.idx.slice      == (slice_nr - 1))      & ...
                        (raw_data.head.idx.contrast   == (contrast_nr - 1))   & ...
                        (raw_data.head.idx.phase      == (phase_nr - 1))      & ...
                        (raw_data.head.idx.repetition == (repetition_nr - 1)) & ...
                        (raw_data.head.idx.set        == (set_nr - 1))        & ...
                        (raw_data.head.idx.segment    == (segment_nr - 1)));

    %% Get k-space data (Nk x Ni x Nc)
    kspace = complex(zeros(Nk, Ni, Nc, 'single'));
    for idx1 = 1:length(profile_list)
        interleaf_nr = raw_data.head.idx.kspace_encode_step_1(profile_list(idx1)) + 1;
        kspace(:,interleaf_nr,:) = kspace(:,interleaf_nr,:) + reshape(raw_data.data{profile_list(idx1)}(index_range,:), [Nk 1 Nc]);
    end

    %% Prewhiten k-space data
    tic; fprintf('Prewhitening k-space data... ');
    kspace = ipermute(reshape(inv_L * reshape(permute(kspace, [3 1 2]), [Nc Nk*Ni]), [Nc Nk Ni]), [3 1 2]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

    %% Calculate the actual slice number for Siemens interleaved multislice imaging
    if nr_slices > 1 % multi-slice
        if mod(nr_slices,2) == 0 % even
            offset1 = 0;
            offset2 = 1;
        else % odd
            offset1 = 1;
            offset2 = 0;
        end
        if slice_nr <= ceil(nr_slices / 2)
            actual_slice_nr = 2 * slice_nr - offset1;
        else
            actual_slice_nr = 2 * (slice_nr - ceil(nr_slices/2)) - offset2;
        end
    else
        actual_slice_nr = slice_nr;
    end

    %% Get a slice normal vector from Siemens raw data
    %----------------------------------------------------------------------
    % dNormalSag: Sagittal component of a slice normal vector (in PCS)
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal, 'dSag')
        dNormalSag = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal.dSag;
    else
        dNormalSag = 0;
    end

    %----------------------------------------------------------------------
    % dNormalCor: Coronal component of a slice normal vector (in PCS)
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal, 'dCor')
        dNormalCor = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal.dCor;
    else
        dNormalCor = 0;
    end

    %----------------------------------------------------------------------
    % dNormalTra: Transverse component of a slice normal vector (in PCS)
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal, 'dTra')
        dNormalTra = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sNormal.dTra;
    else
        dNormalTra = 0;
    end

    %----------------------------------------------------------------------
    % dRotAngle: Slice rotation angle ("swap Fre/Pha")
    %----------------------------------------------------------------------
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}, 'dInPlaneRot')
        dRotAngle = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.dInPlaneRot; % [rad]
    else
        dRotAngle = 0; % [rad]
    end

    %% Get a slice offset in PCS from Siemens raw data
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}, 'sPosition')
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dSag')
            sag_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dSag; % [mm]
        else
            sag_offset = 0; % [mm]
        end
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dCor')
            cor_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dCor; % [mm]
        else
            cor_offset = 0; % [mm]
        end
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition, 'dTra')
            tra_offset = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_nr}.sPosition.dTra; % [mm]
        else
            tra_offset = 0; % [mm]
        end
    else
        sag_offset = 0; % [mm]
        cor_offset = 0; % [mm]
        tra_offset = 0; % [mm]
    end
    PCS_offset = [sag_offset; cor_offset; tra_offset] * 1e-3; % [mm] * [m/1e3mm] => [m]

    main_orientation = fGSLClassOri(dNormalSag, dNormalCor, dNormalTra);
    if main_orientation == SAGITTAL
        PCS_offset(3) = 0;
    end

    %% Calculate spatial coordinates in DCS [m]
    %----------------------------------------------------------------------
    % Calculate a transformation matrix from RCS to GCS [r,c,s] <=> [PE,RO,SL] 
    %----------------------------------------------------------------------
    rotMatrixRCSToGCS = [0    1    0 ; % [PE]   [0 1 0] * [r]
                         1    0    0 ; % [RO] = [1 0 0] * [c]
                         0    0    1]; % [SL]   [0 0 1] * [s]

    %----------------------------------------------------------------------
    % Calculate a rotation matrix from GCS to PCS from Siemens raw data
    %----------------------------------------------------------------------
    [rotMatrixGCSToPCS,PE_sign,RO_sign,main_orientation] = calcMatrixGCSToPCS(dNormalSag, dNormalCor, dNormalTra, dRotAngle);

    %----------------------------------------------------------------------
    % Calculate a rotation matrix from PCS to DCS
    %----------------------------------------------------------------------
    rotMatrixPCSToDCS = calcMatrixPCSToDCS(patient_position);

    %----------------------------------------------------------------------
    % Calculate a rotation matrix from GCS to DCS
    %----------------------------------------------------------------------
    rotMatrixGCSToDCS = rotMatrixPCSToDCS * rotMatrixGCSToPCS;

    %----------------------------------------------------------------------
    % Calculate a scaling matrix
    %----------------------------------------------------------------------
    scaling_matrix = diag(field_of_view_mm ./ matrix_size * 1e-3); % [mm] * [m/1e3mm] => [m]

    %----------------------------------------------------------------------
    % Calculate a slice offset in DCS
    %----------------------------------------------------------------------
    DCS_offset = rotMatrixPCSToDCS * PCS_offset; % 3 x 1

    %----------------------------------------------------------------------
    % Calculate spatial coordinates in DCS
    %----------------------------------------------------------------------
    r_range = (-floor(N1/2):ceil(N1/2)-1).';
    c_range = (-floor(N2/2):ceil(N2/2)-1).';
    s_range = (-floor(N3/2):ceil(N3/2)-1).';
    [r,c,s] = ndgrid(r_range, c_range, s_range);

    xyz = rotMatrixPCSToDCS * rotMatrixGCSToPCS * rotMatrixRCSToGCS * scaling_matrix * cat(2, r(:), c(:), s(:)).';
    x = reshape(xyz(1,:), [N 1]) + DCS_offset(1); % N x 1 [m]
    y = reshape(xyz(2,:), [N 1]) + DCS_offset(2); % N x 1 [m]
    z = reshape(xyz(3,:), [N 1]) + DCS_offset(3); % N x 1 [m]

    %% Display slice information
    fprintf('======================= SLICE INFORMATION =======================\n');
    fprintf('main_orientation = %d (SAGITTAL/CORONAL/TRANSVERSE = 0/1/2)\n', main_orientation);
    fprintf('slice_nr = %d, actual_slice_nr = %d\n', slice_nr, actual_slice_nr);
    fprintf('dNormalSag = %+g \ndNormalCor = %+g \ndNormalTra = %+g \ndRotAngle = %g [rad]\n', dNormalSag, dNormalCor, dNormalTra, dRotAngle);
    fprintf('---------------------- From Siemens raw data --------------------\n');
    fprintf('                   [Sag]   %10.5f [mm]\n', sag_offset);
    fprintf('slice offset(PCS): [Cor] = %10.5f [mm]\n', cor_offset);
    fprintf('                   [Tra]   %10.5f [mm]\n', tra_offset);
    fprintf('---------------------- From Siemens raw data --------------------\n');
    fprintf('                   [Sag]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS(1,1), rotMatrixGCSToPCS(1,2), rotMatrixGCSToPCS(1,3));
    fprintf('rotMatrixGCSToPCS: [Cor] = [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS(2,1), rotMatrixGCSToPCS(2,2), rotMatrixGCSToPCS(2,3));
    fprintf('                   [Tra]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToPCS(3,1), rotMatrixGCSToPCS(3,2), rotMatrixGCSToPCS(3,3));
    fprintf('-----------------------------------------------------------------\n');
    fprintf('                   [Sag]   [%10.5f %10.5f %10.5f]\n', rotMatrixPCSToDCS(1,1), rotMatrixPCSToDCS(1,2), rotMatrixPCSToDCS(1,3));
    fprintf('rotMatrixPCSToDCS: [Cor] = [%10.5f %10.5f %10.5f]\n', rotMatrixPCSToDCS(2,1), rotMatrixPCSToDCS(2,2), rotMatrixPCSToDCS(2,3));
    fprintf('                   [Tra]   [%10.5f %10.5f %10.5f]\n', rotMatrixPCSToDCS(3,1), rotMatrixPCSToDCS(3,2), rotMatrixPCSToDCS(3,3));
    fprintf('-----------------------------------------------------------------\n');
    fprintf('                   [Sag]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToDCS(1,1), rotMatrixGCSToDCS(1,2), rotMatrixGCSToDCS(1,3));
    fprintf('rotMatrixGCSToDCS: [Cor] = [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToDCS(2,1), rotMatrixGCSToDCS(2,2), rotMatrixGCSToDCS(2,3));
    fprintf('                   [Tra]   [%10.5f %10.5f %10.5f]\n', rotMatrixGCSToDCS(3,1), rotMatrixGCSToDCS(3,2), rotMatrixGCSToDCS(3,3));
    fprintf('-----------------------------------------------------------------\n');
    fprintf('                   [ x ]   %10.5f [mm]\n', DCS_offset(1) * 1e3);
    fprintf('slice offset(DCS): [ y ] = %10.5f [mm]\n', DCS_offset(2) * 1e3);
    fprintf('                   [ z ]   %10.5f [mm]\n', DCS_offset(3) * 1e3);
    fprintf('=================================================================\n');

    %% Calculate nominal gradient waveforms in GCS [mT/m]
    %----------------------------------------------------------------------
    % Calculate one spiral interleaf: k in [cycle/cm], g in [G/cm]
    %----------------------------------------------------------------------
    Fcoeff = [FOV -FOV * (1 - user_opts.vds_factor / 100)]; % FOV decreases linearly from 30 to 15cm
    krmax  = 1 / (2 * (FOV / matrix_size(1)));              % [cycle/cm]
    res    = 1 / (2 * krmax) * 10;                          % resolution [mm]
    [k_spiral_arm,g_spiral_arm,s_spiral_arm,time] = vdsmex(interleaves, Fcoeff, res, Gmax, Smax, sampling_time, 10000000);
    g_spiral_arm = -g_spiral_arm(1:Nk,:); % Nk x 2

    %----------------------------------------------------------------------
    % Rotate the spiral interleaf by 360/Ni degrees every TR
    % (Re{g_spiral} + 1j * Im{g_spiral}) * (cos(arg) - 1j * sin(arg))
    % RO: real =>  Re{g_spiral} * cos(arg) + Im{g_spiral} * sin(arg)
    % PE: imag => -Re{g_spiral} * sin(arg) + Im{g_spiral} * cos(arg)
    %----------------------------------------------------------------------
    gu_nominal = zeros(Nk, Ni, 'double'); % Nk x Ni [G/cm] PE
    gv_nominal = zeros(Nk, Ni, 'double'); % Nk x Ni [G/cm] RO
    gw_nominal = zeros(Nk, Ni, 'double'); % Nk x Ni [G/cm] SL
    for i = 1:Ni
        arg = 2 * pi / Ni * (i - 1); % [rad]
        gv_nominal(:,i) =  g_spiral_arm(:,1) * cos(arg) + g_spiral_arm(:,2) * sin(arg);
        gu_nominal(:,i) = -g_spiral_arm(:,1) * sin(arg) + g_spiral_arm(:,2) * cos(arg);
    end

    %% Calculate GIRF-predicted gradient waveforms [rad/m]
    g_nominal = cat(3, gu_nominal, gv_nominal);
    tRR = 0; % custom clock-shift
    sR.R = rotMatrixGCSToDCS;
    sR.T = header.acquisitionSystemInformation.systemFieldStrength_T;
    [~,g_predicted] = apply_GIRF(g_nominal, dt, sR, tRR); % k:[cycle/cm] and g:[G/cm]

    %% Change the sign of GIRF-predicted gradient waveforms [G/cm]
    gu_predicted = PE_sign * g_predicted(:,:,1); % [G/cm] PE
    gv_predicted = RO_sign * g_predicted(:,:,2); % [G/cm] RO
    gw_predicted = g_predicted(:,:,3);           % [G/cm] SL

    %% Calculate GIRF-predicted k-space trajectories in GCS [rad/m]
    %----------------------------------------------------------------------
    % Numerically integrate the coefficients
    % [rad/sec/T] * [G/cm] * [T/1e4G] * [1e2cm/m] * [sec] => [rad/m]
    %----------------------------------------------------------------------
    ku_predicted = cumsum(gamma * gu_predicted * 1e-2 * double(dt)); % [rad/m] PE
    kv_predicted = cumsum(gamma * gv_predicted * 1e-2 * double(dt)); % [rad/m] RO
    kw_predicted = cumsum(gamma * gw_predicted * 1e-2 * double(dt)); % [rad/m] SL

    %% Transform gradient waveforms from GCS to DCS
    gx_predicted = zeros(Nk, Ni, 'double'); % [G/cm]
    gy_predicted = zeros(Nk, Ni, 'double'); % [G/cm]
    gz_predicted = zeros(Nk, Ni, 'double'); % [G/cm]
    for i = 1:Ni
        g_log = cat(2, gu_predicted(:,i), gv_predicted(:,i), gw_predicted(:,i)).'; % 3 x Nk
        g_phy = rotMatrixGCSToDCS * g_log; % 3 x Nk
        gx_predicted(:,i) = g_phy(1,:); % Nk x 1
        gy_predicted(:,i) = g_phy(2,:); % Nk x 1
        gz_predicted(:,i) = g_phy(3,:); % Nk x 1
    end

    %% Calculate concomitant field basis functions (N x Nl) [m], [m^2], [m^3]
    tic; fprintf('Calculating concomitant field basis functions... ');
    p = calculate_concomitant_field_basis(x, y, z, Nl);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

    %% Calculate the time courses of phase coefficients (Nk x Nl x Ni) [rad/m], [rad/m^2], [rad/m^3]
    tic; fprintf('Calculating the time courses of phase coefficients... ');
    k = calculate_concomitant_field_coefficients(gx_predicted, gy_predicted, gz_predicted, Nl, B0, gamma, dt);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    
    %% Calculate a density compensation function
    local_tic = tic; fprintf('Cacluating a density compensation function based on sdc3_MAT.c... ');
    % [rad/m] / ([cycle/cm] * [2pi rad/cycle] * [1e2cm/m]) => [unitless]
    coords = cat(1, reshape(kv_predicted   / (2 * krmax * 2 * pi * 1e2), [1 Nk Ni]), ...
                    reshape(ku_predicted   / (2 * krmax * 2 * pi * 1e2), [1 Nk Ni]), ...
                    reshape(0*kw_predicted / (2 * krmax * 2 * pi * 1e2), [1 Nk Ni]));
    numIter = 25;
    effMtx  = matrix_size(1);
    verbose = 0;
    DCF = sdc3_MAT(coords, numIter, effMtx, verbose, 2.1);
    w = DCF / max(DCF(:));
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(local_tic), toc(start_time));

    %% Calculate a circular mask (image support constraint)
    switch main_orientation % (S/C/T = 0/1/2)
        case SAGITTAL
            mask = reshape((y.^2 + z.^2 < (field_of_view_mm(1) * osf * 1e-3 / 2)^2), [N1 N2]); % [mm] * [m/1e3mm] => [m]
        case CORONAL
            mask = reshape((x.^2 + z.^2 < (field_of_view_mm(1) * osf * 1e-3 / 2)^2), [N1 N2]); % [mm] * [m/1e3mm] => [m]
        case TRANSVERSE
            mask = reshape((x.^2 + y.^2 < (field_of_view_mm(1) * osf * 1e-3 / 2)^2), [N1 N2]); % [mm] * [m/1e3mm] => [m]
        otherwise
    end

    %----------------------------------------------------------------------
    % Combine a circular mask with a support of B0 maps
    %----------------------------------------------------------------------
    %B0mask = (abs(B0map) > 0);
    %mask = mask | B0mask;

    %% Perform NUFFT reconstruction
    %----------------------------------------------------------------------
    % Prepare an NUFFT structure using k-space trajectories in RCS
    % Note: RCS <=> GCS (PRS)
    % ku:[PE]   [0 1 0] * [r]
    % kv:[RO] = [1 0 0] * [c]
    % kw:[SL]   [0 0 1] * [s]
    %----------------------------------------------------------------------
    % scaled to [-0.5,0.5] and then [-pi,pi]
    % [rad/m] / ([cycle/cm] * [2pi rad/cycle] * [1e2cm/m]) => [rad/m] / [rad/m] = [unitless]
    k_rcs = cat(2, reshape(kv_predicted, [Nk*Ni 1]), reshape(ku_predicted, [Nk*Ni 1])) / (2 * krmax * 2 * pi * 1e2); % [-0.5,0.5]
    om = -k_rcs * 2 * pi; % The definition of FFT is opposite in NUFFT
    Nd = [N1 N2]; % matrix size
    Jd = [6 6];   % kernel size
    Kd = Nd * 2;  % oversampled matrix size
    nufft_st = nufft_init(om, Nd, Jd, Kd, Nd/2, 'minmax:kb');

    start_time_nufft = tic;
    imc_nufft = complex(zeros(N1, N2, Nc, 'double'));
    for c = 1:Nc
        tic; fprintf('(%2d/%2d): NUFFT reconstruction (c=%2d/%2d)... ', idx2, nr_recons, c, Nc);
        %------------------------------------------------------------------
        % Apply the adjoint of 2D NUFFT
        %------------------------------------------------------------------
        imc_nufft(:,:,c) = nufft_adj(reshape(double(kspace(:,:,c)) .* w, [Nk*Ni 1]), nufft_st) / sqrt(prod(Nd));
        fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    end
    computation_time_nufft = toc(start_time_nufft);

    %% Calculate coil sensitivity maps
    tic; fprintf('(%2d/%2d): Calculating coil sensitivity maps with Walsh method... ', idx2, nr_recons);
    %----------------------------------------------------------------------
    % IFFT to k-space (k-space <=> image-space)
    %----------------------------------------------------------------------
    kspace_gridded = ifft2c(imc_nufft);

    %----------------------------------------------------------------------
    % Calculate the calibration region of k-space
    %----------------------------------------------------------------------
    cal_shape = [32 32];
    cal_data = crop(kspace_gridded, [cal_shape Nc]);
    cal_data = bsxfun(@times, cal_data, hamming(cal_shape(1)) * hamming(cal_shape(2)).');

    %----------------------------------------------------------------------
    % Calculate coil sensitivity maps
    %----------------------------------------------------------------------
    cal_im = fft2c(zpad(cal_data, [N1 N2 Nc]));
    csm = ismrm_estimate_csm_walsh(cal_im);

    %----------------------------------------------------------------------
    % Apply a circular mask to CSM
    %----------------------------------------------------------------------
    csm = bsxfun(@times, csm, mask);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

    %% Perform optimal coil combination for NUFFT reconstruction
    im_nufft = sum(imc_nufft .* conj(csm), 3);
    figure, imagesc(abs(im_nufft).'); axis image; colormap(gray(256));

    %% Perform CG-SENSE reconstruction
    start_time_sense = tic;
    tic; fprintf('(%2d/%2d): Performing CG-SENSE reconstruction...\n', idx2, nr_recons);
    b = reshape(sum(conj(csm).* imc_nufft, 3), [N 1]);
    max_iterations = 30;
    limit = 1e-5;
    E = @(x,tr) encoding_SENSE(x, csm, w, nufft_st, tr);
    [m_sense, flag, relres, iter, resvec, lsvec] = lsqr(E, b, limit, max_iterations, [], [], []); % N x 1
    im_sense = reshape(m_sense, [N1 N2]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    computation_time_sense = toc(start_time_sense);

    %% Estimate coil sensitivity maps using ESPIRiT
    if support_constraint
        tic; fprintf('(%2d/%2d): Calculating coil sensitivity maps with ESPIRiT... ', idx2, nr_recons);
        %------------------------------------------------------------------
        % Set parameters
        %------------------------------------------------------------------
        ncalib = 24;    % use 24 calibration lines to compute compression
        ksize  = [6 6]; % kernel size

        % Threshold for picking singular vercors of the calibration matrix
        % (relative to largest singlular value.
        eigThresh_1 = 0.02;

        % threshold of eigen vector decomposition in image space.
        eigThresh_2 = 0.95;

        %------------------------------------------------------------------
        % Calculate "clean" k-space data
        %------------------------------------------------------------------
        kspace_sense = fft2c(bsxfun(@times, csm, im_sense));

        %------------------------------------------------------------------
        % crop a calibration area
        %------------------------------------------------------------------
        calib = crop(kspace_sense, [ncalib ncalib Nc]);

        %------------------------------------------------------------------
        % Compute ESPIRiT EigenVectors
        % We perform calibration in k-space followed by an eigen-decomposition
        % in image space to produce the EigenMaps.
        %------------------------------------------------------------------
        % compute the calibration matrix, perform SVD, and convert singular
        % vectors into k-space kernels
        [kernel,S] = dat2Kernel(calib, ksize);
        index = max(find(S >= S(1) * eigThresh_1));

        %------------------------------------------------------------------
        % Compute eigen-value decomposition in image space to get CSM maps
        %------------------------------------------------------------------
        [V,W] = kernelEig(kernel(:,:,:,1:index), [N1 N2]);

        %------------------------------------------------------------------
        % Calculate a "support" mask based on an eigenvalue threshold
        %------------------------------------------------------------------
        mask_espirit = (W(:,:,end) > eigThresh_2);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));

        %% Apply a support mask to CSM
        csm = bsxfun(@times, csm, mask_espirit);
    end

    %% Perform King's method for concomitant field correction
    start_time_king = tic;
    tic; fprintf('(%2d/%2d): Performing King''s method...\n', idx2, nr_recons);
    [im_king,im_fs] = perform_deblurring_king_method(kspace, nufft_st, w, csm, gx_predicted, gy_predicted, gz_predicted, x, y, z, rotMatrixRCSToGCS, rotMatrixGCSToPCS, rotMatrixPCSToDCS, field_of_view_mm, DCS_offset, gamma, B0, dt);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    computation_time_king = toc(start_time_king);

    %% Calculate a time vector [sec]
    t = TE + (0:Nk-1).' * dt; % Nk x 1 [sec]

    %% Calculate the SVD of the higher-order encoding matrix (Nk x N)
    os = 5; % oversampling parameter for randomized SVD
    u_tilde = complex(zeros(Nk, Lmax, Ni, 'double'));
    v_tilde = complex(zeros(N , Lmax, Ni, 'double'));
    singular_values = complex(zeros(Lmax, Lmax, Ni, 'double'));
    for i = 1:Ni
        tstart = tic; fprintf('(%2d/%2d): Calculating the randomized SVD (i=%2d/%2d)... \n', idx2, nr_recons, i, Ni);
        [U,S,V] = calculate_rsvd_higher_order_encoding_matrix(k(:,4:end,i), p(:,4:end), Lmax, os, vec(B0map(:,:,actual_slice_nr)), t, static_B0_correction);
        u_tilde(:,:,i) = U(:,1:Lmax); % Nk x Lmax
        v_tilde(:,:,i) = V(:,1:Lmax) * S(1:Lmax,1:Lmax)'; % N x Lmax
        singular_values(:,:,i) = S(1:Lmax,1:Lmax);
        fprintf('(%2d/%2d): Calculating the randomized SVD (i=%2d/%2d)... done! (%6.4f/%6.4f sec)\n', idx2, nr_recons, i, Ni, toc(tstart), toc(start_time));
    end

    %% Calculate NUFFT structures for lowrank MaxGIRF reconstruction
    %----------------------------------------------------------------------
    % Prepare an NUFFT structure using k-space trajectories in RCS
    % Note: RCS <=> GCS (PRS)
    % ku:[PE]   [0 1 0] * [r]
    % kv:[RO] = [1 0 0] * [c]
    % kw:[SL]   [0 0 1] * [s]
    %----------------------------------------------------------------------
    st = cell(Ni,1);
    for i = 1:Ni
        tic; fprintf('(%2d/%2d): Calculating NUFFT structure... ', i, Ni);
        % scaled to [-0.5,0.5] and then [-pi,pi]
        % [rad/m] / ([cycle/cm] * [2pi rad/cycle] * [1e2cm/m]) => [rad/m] / [rad/m] = [unitless]
        k_rcs = cat(2, kv_predicted(:,i), ku_predicted(:,i)) / (2 * krmax * 2 * pi * 1e2); % [-0.5,0.5]
        om = -k_rcs * 2 * pi; % The definition of FFT is opposite in NUFFT
        Nd = [N1 N2];   % matrix size
        Jd = [6 6];     % kernel size
        Kd = Nd * 2;    % oversampled matrix size
        st{i} = nufft_init(om, Nd, Jd, Kd, Nd/2, 'minmax:kb');
        fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    end

    %% Perform conjugate phase reconstruction using the MaxGIRF encoding model
    b = zeros(N, Lmax, 'double');
    for i = 1:Ni
        Nd = st{i}.Nd;
        tic; fprintf('(%d/%d): Calculating conjugate phase reconstruction (i=%d/%d)... ', idx2, nr_recons, i, Ni);

        for c = 1:Nc
            %--------------------------------------------------------------
            % Caclulate d_{i,c}
            %--------------------------------------------------------------
            d = double(kspace(:,i,c)); % kspace: Nk x Ni x Nc

            %--------------------------------------------------------------
            % Calculate sum_{ell=1}^L ...
            % diag(v_tilde(ell,i)) * Fc^H * diag(conj(u_tilde(ell,i)))
            %--------------------------------------------------------------
            AHd  = zeros(N, Lmax, 'double');
            AHd_ = zeros(N, 1, 'double');
            for ell = 1:Lmax
                % Preconditioning with density compensation
                FHDuHd = nufft_adj((conj(u_tilde(:,ell,i)) .* d) .* w(:,i), st{i}) / sqrt(prod(Nd));
                AHd_ = AHd_ + v_tilde(:,ell,i) .* reshape(FHDuHd, [N 1]);
                AHd(:,ell) = AHd_;
            end

            %--------------------------------------------------------------
            % Calculate Sc^H * Ei^H * d_{i,c}
            %--------------------------------------------------------------
            for ell = 1:Lmax
                AHd(:,ell) = reshape(conj(double(csm(:,:,c))), [N 1]) .* AHd(:,ell);
            end

            %--------------------------------------------------------------
            % Calculate b (N x 1)
            %--------------------------------------------------------------
            for ell = 1:Lmax
                b(:,ell) = b(:,ell) + AHd(:,ell);
            end
        end
        fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    end
    im_cpr = reshape(b, [N1 N2 Lmax]);

    %% Perform lowrank MaxGIRF reconstruction
    start_time_lowrank = tic;
    tstart = tic; fprintf('(%d/%d): Performing lowrank MaxGIRF reconstruction...\n', idx2, nr_recons);
    max_iterations = 15;
    limit = 1e-5;
    E = @(x,tr) encoding_lowrank_MaxGIRF(x, csm, u_tilde(:,1:L,:), v_tilde(:,1:L,:), w, st, tr);
    [m_lowrank, flag, relres, iter, resvec, lsvec] = lsqr(E, b(:,L), limit, max_iterations, [], [], []); % NL x 1
    im_lowrank = reshape(m_lowrank, [N1 N2]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    computation_time_lowrank = toc(start_time_lowrank);

    return
    %% Case 1: MaxGIRF (predicted, without off-resonance, without conc. field correction) 
    E = complex(zeros(Nk, N, Ni, 'single'));
    for i = 1:Ni
        tic; fprintf('(i=%d/%d): Calculating the encoding matrix E... ', i, Ni);
        % (Nk x 3) x (N x 3).' => Nk x N
        phi = k(:,1:3,i) * p(:,1:3).';
        E(:,:,i) = exp(1j * phi);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    end

    %----------------------------------------------------------------------
    % Iterative image reconstruction
    %----------------------------------------------------------------------
    start_time_maxgirf1 = tic;
    tic; fprintf('Performing MaxGIRF reconstruction...\n');
    [im1,delta1] = iterative_image_reconstruction_precompute(kspace, E, single(csm), beta, maxiter);
    im_maxgirf1 = reshape(im1(:,maxiter), [N1 N2]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    computation_time_maxgirf1 = toc(start_time_maxgirf1);

    %% Case 2: MaxGIRF (predicted, without off-resonance, with conc. field correction) 
    E = complex(zeros(Nk, N, Ni, 'single'));
    for i = 1:Ni
        tic; fprintf('(i=%d/%d): Calculating the encoding matrix E... ', i, Ni);
        % (Nk x 3) x (N x 3).' => Nk x N
        phi = k(:,:,i) * p.';
        E(:,:,i) = exp(1j * phi);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    end

    %----------------------------------------------------------------------
    % Iterative image reconstruction
    %----------------------------------------------------------------------
    start_time_maxgirf2 = tic;
    tic; fprintf('Performing MaxGIRF reconstruction...\n');
    [im2,delta2] = iterative_image_reconstruction_precompute(kspace, E, single(csm), beta, maxiter);
    im_maxgirf2 = reshape(im2(:,maxiter), [N1 N2]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    computation_time_maxgirf2 = toc(start_time_maxgirf2);

    %% Case 3: MaxGIRF (predicted, with off-resonance, with conc. field correction) 
    E = complex(zeros(Nk, N, Ni, 'single'));
    for i = 1:Ni
        tic; fprintf('(i=%d/%d): Calculating the encoding matrix E... ', i, Ni);
        % (Nk x Nl) x (N x Nl).' => Nk x N
        phi = 2 * pi * t * B0map(:).' + k(:,:,i) * p.';
        E(:,:,i) = exp(1j * phi);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    end

    %----------------------------------------------------------------------
    % Iterative image reconstruction
    %----------------------------------------------------------------------
    start_time_maxgirf3 = tic;
    tic; fprintf('Performing MaxGIRF reconstruction...\n');
    [im3,delta3] = iterative_image_reconstruction_precompute(kspace, E, single(csm), beta, maxiter);
    im_maxgirf3 = reshape(im3(:,maxiter), [N1 N2]);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    computation_time_maxgirf3 = toc(start_time_maxgirf3);
end


%% Calculate a time-varying concomitant field map
KP = complex(zeros(Nk, N, Ni, 'single'));
for i = 1:Ni
    tic; fprintf('(i=%d/%d): Calculating the KP... ', i, Ni);
    % (Nk x 3) x (N x 3).' => Nk x N
    KP(:,:,i) = k(:,4:end,i) * p(:,4:end).';
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end

if 0
%% Make a movie
X = reshape(x*1e3, [N1 N2]);
Y = reshape(y*1e3, [N1 N2]);
Z = reshape(z*1e3, [N1 N2]);

figure('Color', 'k', 'Position', [3 2 425 770]);

% create video writer object
vidfile = VideoWriter('spiral_and_concomitant_field_map.mp4', 'MPEG-4');

% set the frame rate to one frame per second
set(vidfile, 'FrameRate', 5);
open(vidfile);

text_color = 'w';
color_order = get(gca, 'colororder');

ax1 = subplot(3,1,1); hold on;
ax2 = subplot(3,1,2); hold on;
ax3 = subplot(3,1,3); hold on;

pos2 = get(ax2, 'Position');
pos3 = get(ax3, 'Position');
set(ax1, 'Position', [0.1500 0.7380 0.7494 0.1683]); % [left bottom width height]
set(ax2, 'Position', [0.2300 0.3300 0.6293 0.2500]); % [left bottom width height]
set(ax3, 'Position', [0.2300 0.0200 0.6291 0.2500]); % [left bottom width height]
hText = text(ax2, 0, 0, sprintf('Concomitant field map at t = %3.2f msec', t(idx1)*1e3), 'Color', text_color, 'rotation', 0, 'FontSize', 14);
set(hText, 'Position', [-53.0570 198.9637 0], 'FontWeight', 'Bold');

for idx2 = 1:Ni
    hPlot1 = plot(ax1, t*1e3, gx_predicted(:,idx2)*1e1, 'LineWidth', 2, 'Color', color_order(1,:));
    hPlot2 = plot(ax1, t*1e3, gy_predicted(:,idx2)*1e1, 'LineWidth', 2, 'Color', color_order(2,:));
    h_current(1) = line(ax1, t(1)*1e3, abs(gx_predicted(1,idx2)*1e1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
    h_current(2) = line(ax1, t(1)*1e3, abs(gy_predicted(1,idx2)*1e1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'r');
    h_current(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
    h_current(2).Annotation.LegendInformation.IconDisplayStyle = 'off';

    xlabel(ax1, 'Time [msec]', 'Color', text_color);
    ylabel(ax1, 'Amplitude [mT/m]', 'Color', text_color);
    %hLegend = legend(ax1, '$G_x$', '$G_y$', 'Interpreter', 'latex', 'Color', 'k', 'TextColor', text_color, 'FontSize', 14);
    %set(hLegend, 'Orientation', 'horizontal', 'Position', [0.0902 0.6114 0.1743 0.0640]);
    title(ax1, {'GIRF-predicted phyiscal gradients', sprintf('Interleaf number = %d', idx2)}, 'Color', text_color);
    xlim(ax1, [0 t(end)]*1e3);
    set(ax1, 'FontSize', 14, 'Color', 'k', 'XColor', text_color, 'YColor', text_color);

    for idx1 = 1:80:Nk
        % Update graphics data. This is more efficient than recreating plots.
        set(hText, 'string', sprintf('Concomitant field map at t = %3.2f msec', t(idx1)*1e3));
        set(h_current(1), 'XData', t(idx1)*1e3, 'YData', gx_predicted(idx1,idx2)*1e1);
        set(h_current(2), 'XData', t(idx1)*1e3, 'YData', gy_predicted(idx1,idx2)*1e1);

        imagesc(ax2, reshape(KP(idx1,:,idx2)*0, [N1 N2]));
        axis(ax2, 'image', 'off');
        title(ax2, sprintf('Isocenter'), 'Color', text_color);
        set(ax2, 'FontSize', 14, 'XColor', text_color, 'YColor', text_color);
        colormap(jet(256));
        hc = colorbar(ax2);
        caxis(ax2, [0 100]);
        set(hc, 'Color', text_color);
        title(hc, '[Hz]', 'Color', text_color);

        imagesc(ax3, reshape(KP(idx1,:,idx2) / (2 * pi * t(end)), [N1 N2]));
        axis(ax3, 'image', 'off');
        title(ax3, sprintf('z = 75 mm'), 'Color', text_color);
        set(ax3, 'FontSize', 14, 'XColor', text_color, 'YColor', text_color);
        colormap(ax3, jet(256));
        hc = colorbar(ax3);
        caxis(ax3, [0 100]);
        set(hc, 'Color', text_color);
        title(hc, '[Hz]', 'Color', text_color);
        drawnow;

        frame = getframe(gcf);
        writeVideo(vidfile, frame);
    end
    delete(hPlot1);
    delete(hPlot2);    
    delete(h_current(1));
    delete(h_current(2));
end
close(vidfile);
end
    %%
    
    
    im1 = flip(rot90(abs(im_nufft) / max(abs(im_nufft(:))),-1),2);
    im2 = flip(rot90(abs(im_sense) / max(abs(im_sense(:))),-1),2);
    im3 = flip(rot90(abs(im_king) / max(abs(im_king(:))),-1),2);
    im4 = flip(rot90(abs(im_lowrank) / max(abs(im_lowrank(:))),-1),2);
    im5 = flip(rot90(abs(im_maxgirf1) / max(abs(im_maxgirf1(:))),-1),2);
    im6 = flip(rot90(abs(im_maxgirf2) / max(abs(im_maxgirf2(:))),-1),2);
    im7 = flip(rot90(abs(im_maxgirf3) / max(abs(im_maxgirf3(:))),-1),2);
    im8 = flip(rot90(abs(im_gre(:,:,1,1)) / max(abs(vec(im_gre(:,:,1,1)))),-1),2);   
    im_montage1 = cat(2, im1  , im2, im3, im4, im5*0);
    im_montage2 = cat(2, im1*0, im5, im6, im7, im8);
    im_montage = cat(1, im_montage1, im_montage2);
    figure('Color', 'w');
    imagesc(abs(im_montage)); axis image; colormap(gray(256));
    caxis([0 0.7]);


    
    
    

return
%%

figure('Color', 'w');
for idx2 = 1:maxiter
    imagesc(flip(rot90(abs(reshape(im3(:,idx2), [N1 N2])),-1),2)); axis image;
    colormap(gray(256));
    colorbar;
    title(sprintf('Iteration number = %d', idx2));
    pause;
end
%%
figure('Color', 'w');
plot((1:maxiter).', delta, '.-', 'MarkerSize', 10);
set(gca, 'YScale', 'log');

%%
figure; montage(flip(rot90(reshape(abs(im_maxgirf3), [N1 N2 maxiter]),-1),2), 'DisplayRange', []);







%% Display results
mask = (abs(B0map) > 0);
im_reference = flip(rot90(im_gre(:,:,1),-1),2) / max(abs(vec(im_gre(:,:,1))));
im1 = flip(rot90(im_maxgirf1,-1),2) / max(abs(vec(im_maxgirf1)));
im2 = flip(rot90(im_maxgirf2,-1),2) / max(abs(vec(im_maxgirf2)));
im3 = flip(rot90(im_maxgirf3,-1),2) / max(abs(vec(im_maxgirf3)));
im_montage = cat(2, im1, im2, im3, im_reference);

figure('Color', 'w', 'Position', [1 184 1598 615]);
imagesc(abs(im_montage)); axis image off;
colormap(gray(256));
colorbar;

figure('Color', 'w', 'Position', [1 184 1598 615]);
imagesc(angle(im_montage)); axis image off;
colormap(hsv(256));
colorbar;
%%
figure('Color', 'w');
imagesc(flip(rot90(B0map,-1),2)); axis image off;
colormap(jet(256));
colorbar;


return
%%

% %%
% figure;
% for idx = 1:maxiter
%     imagesc(abs(im_maxgirf(:,:,idx))); axis image; colormap(gray(256));
%     title(sprintf('maxiter = %d', idx));
%     pause;
% end

%     maxiter = 20;
%     im_maxgirf = complex(zeros(N1, N2, maxiter, 'single'));
%     for idx1 = 1:maxiter
%         fprintf('MaxGIRF reconstruction... ');
%         im = iterative_image_reconstruction_blas(kspace / scale_factor, E, csm, beta, idx1);
%         im = reshape(im, [N1 N2]);
%         im_maxgirf(:,:,idx1) = im;
%         fprintf('done! (%5.3f sec)\n', toc(start_time));
%         figure; imagesc(abs(flip(rot90(im,-1),2))); axis image; colormap(gray(256)); drawnow;
%     end


    %%
    im1 = flip(rot90(abs(im_nufft) / max(abs(im_nufft(:))),-1),2);
    im2 = flip(rot90(abs(im_sense) / max(abs(im_sense(:))),-1),2);
    im3 = flip(rot90(abs(im_king) / max(abs(im_king(:))),-1),2);
    im4 = flip(rot90(abs(im_maxgirf) / max(abs(im_maxgirf(:))),-1),2);
    im5 = flip(rot90(abs(im_gre(:,:,1,1)) / max(abs(vec(im_gre(:,:,1,1)))),-1),2);   
    im_montage = cat(2, im1, im2, im3, im4, im5);
    figure('Color', 'w');
    imagesc(abs(im_montage)); axis image; colormap(gray(256));
    
    
    
    %%
    
    figure; imagesc(flip(rot90(mask,-1),2)); axis image;
    figure; imagesc(flip(rot90(B0map,-1),2)); axis image; colormap(jet(256));
    figure; imagesc(flip(rot90(abs(im_gre(:,:,1,1)),-1),2)); axis image; colormap(jet(256));
    
%     figure, imagesc(abs(im_sense).'); axis image; colormap(gray(256));
%     figure, imagesc(abs(im_nufft).'); axis image; colormap(gray(256));
%     figure, imagesc(abs(im_king).'); axis image; colormap(gray(256));
%     figure, imagesc(abs(im_maxgirf).'); axis image; colormap(gray(256));
    
    
    
    return;
    
    
    
    
    %%
    %----------------------------------------------------------------------
    % Iterative image reconstruction per coil
    %----------------------------------------------------------------------
    maxiter_csm = 15;   % maximum number of CG iterations
    beta_csm    = 1e-6; % Tikhonov regularization parameter

    imc = complex(zeros(N1, N2, Nc, 'single'));
    for c = 1:Nc
        tic; fprintf('Performing MaxGIRF reconstruction per coil (%d/%d)...\n', c, Nc);
        [im,delta] = iterative_image_reconstruction_precompute(kspace(:,:,c), E, complex(ones(N1, N2, 'single')), beta_csm, maxiter_csm);
        imc(:,:,c) = reshape(im(:,maxiter), [N1 N2]);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    end

    %----------------------------------------------------------------------
    % Scale k-space data such that the maximum absolute value of all coil 
    % images is "1"
    %----------------------------------------------------------------------
    scale_factor1 = max(abs(imc(:)));

    %%
    imc2 = complex(zeros(N1, N2, Nc, 'single'));
    for c = 1:Nc
        tic; fprintf('Performing MaxGIRF reconstruction per coil (%d/%d)...\n', c, Nc);
        [im,delta] = iterative_image_reconstruction_memory_efficient(kspace(:,:,c), k_t, h, t, B0map, complex(ones(N1, N2, 'single')), beta_csm, maxiter_csm);
        imc2(:,:,c) = reshape(im(:,maxiter), [N1 N2]);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
    end
    scale_factor2 = max(abs(imc2(:)));
    
    
    
    %% Calculate coil sensitivity maps
    %----------------------------------------------------------------------
    % IFFT to k-space (k-space <=> image-space)
    %----------------------------------------------------------------------
    kspace_cartesian = ifft2c(imc); % N1 x N2 x Nc

    %----------------------------------------------------------------------
    % Calculate the calibration region of k-space
    %----------------------------------------------------------------------
    cal_shape = [32 32];
    cal_data = crop(kspace_cartesian, [cal_shape Nc]);
    cal_data = bsxfun(@times, cal_data, hamming(cal_shape(1)) * hamming(cal_shape(2)).');

    %----------------------------------------------------------------------
    % Estimate coil sensitivity maps (N1 x N2 x Nc)
    %----------------------------------------------------------------------
    cal_im = fft2c(zpad(cal_data, [N1 N2 Nc]));
    csm = single(ismrm_estimate_csm_walsh(cal_im));
